#include "PIC.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

//  Surface tension 
//  Label(i,j) : i,j ∈ [0, nx+1] × [0, ny+1]  (ghost cells)
//  phi/kappa/normalX/normalY(i,j) : i,j ∈ [0,nx-1] × [0,ny-1]
//  phi(i,j) <-> Label(i+1, j+1)
//  interface_u(i,j) : (nx+1)×ny - face btw phi(i-1,j) & phi(i,j)
//  interface_v(i,j) : nx×(ny+1) - face btw phi(i,j-1) & phi(i,j)

static const varType CSF_KERNEL_RADIUS_FACTOR = 2.0;

static inline varType smoothingKernel(varType r, varType h)
{
    const varType q = r / h;
    if (q >= 1.0) return 0.0;
    const varType s = 1.0 - q;
    return s * s * s;
}

void PIC::UpdatePhiFromParticles() const
{
    const int     nx       = fields->nx;
    const int     ny       = fields->ny;
    const varType dx       = fields->dx;
    const varType dy       = fields->dy;
    const varType h_smooth = CSF_KERNEL_RADIUS_FACTOR * dx;
    const int     Rsmooth  = static_cast<int>(std::ceil(h_smooth / dx)) + 1;

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fields->phi->Set(i, j, 0.0);
            fields->kappa->Set(i, j, 0.0);
            fields->normalX->Set(i, j, 0.0);
            fields->normalY->Set(i, j, 0.0);
        }
    }

    //  color function phi(i,j) = sum on p [ kernel(|xc - xp|, h)]
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {

            const varType xc = (i + 0.5) * dx;
            const varType yc = (j + 0.5) * dy;

            const int ci_lo = std::max(0, i - Rsmooth);
            const int ci_hi = std::min(nx - 1, i + Rsmooth);
            const int cj_lo = std::max(0, j - Rsmooth);
            const int cj_hi = std::min(ny - 1, j + Rsmooth);

            varType cval = 0.0;
            for (int cj = cj_lo; cj <= cj_hi; ++cj) {
                for (int ci = ci_lo; ci <= ci_hi; ++ci) {
                    const Particles &cell = (*cloud)(ci, cj);
                    const int np = cell.size();
                    for (int p = 0; p < np; ++p) {
                        const varType r = std::hypot(cell.GetX(p) - xc,
                                                     cell.GetY(p) - yc);
                        cval += smoothingKernel(r, h_smooth);
                    }
                }
            }
            fields->phi->Set(i, j, cval);
        }
    }

    // normalisation of phi in [0,1]
    varType Cmax = 0.0;
    OMP_PRAGMA(omp parallel for collapse(2) reduction(max:Cmax) schedule(static))
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            Cmax = std::max(Cmax, fields->phi->Get(i, j));

    if (Cmax < REAL_EPSILON) return;

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            fields->phi->Set(i, j, fields->phi->Get(i, j) / Cmax);

    //   phi < 0 : FLUIDE,  phi > 0 : AIR
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const varType val = fields->phi->Get(i, j);
            fields->phi->Set(i, j, IS_FLUID(fields->Label(i + 1, j + 1))
                                  ? -std::abs(val)
                                  :  std::abs(val));
        }
    }
}

void PIC::ComputeSurfaceTensionOnFaces() const
{
    const int nx = fields->nx;
    const int ny = fields->ny;
    const varType dx = fields->dx;
    const varType dy = fields->dy;
    const varType gamma = params.gamma;

    // back to 0
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < fields->interface_u->ny; ++j)
        for (int i = 0; i < fields->interface_u->nx; ++i)
            fields->interface_u->Set(i, j, 0.0);

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < fields->interface_v->ny; ++j)
        for (int i = 0; i < fields->interface_v->nx; ++i)
            fields->interface_v->Set(i, j, 0.0);

    // normal n = grad(phi)/|grad(phi)| at cell centers grad(phi)
    // using finite differences. Borders clamped to nearest valid cell
    // i.e. d(phi)/d(n) = 0 at domain borders
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {

            const varType phi_xp = fields->phi->Get(std::min(i + 1, nx - 1), j);
            const varType phi_xm = fields->phi->Get(std::max(i - 1, 0),       j);
            const varType phi_yp = fields->phi->Get(i, std::min(j + 1, ny - 1));
            const varType phi_ym = fields->phi->Get(i, std::max(j - 1, 0)      );

            const varType gx = (phi_xp - phi_xm) / (2.0 * dx);
            const varType gy = (phi_yp - phi_ym) / (2.0 * dy);
            const varType gn = std::hypot(gx, gy);

            if (gn > REAL_EPSILON) {
                fields->normalX->Set(i, j, gx / gn);
                fields->normalY->Set(i, j, gy / gn);
            }
        }
    }

    // curvature kappa = -(d(nx)/d(x) + d(ny)/d(y)) at cell centers
    // normalX(i±1, j) & normalY(i, j±1) clamped as above
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {

            if (!IS_FLUID(fields->Label(i + 1, j + 1))) continue;

            const varType dnx_dx =
                (fields->normalX->Get(std::min(i + 1, nx - 1), j) -
                 fields->normalX->Get(std::max(i - 1, 0),       j)) / (2.0 * dx);

            const varType dny_dy =
                (fields->normalY->Get(i, std::min(j + 1, ny - 1)) -
                 fields->normalY->Get(i, std::max(j - 1, 0)      )) / (2.0 * dy);

            fields->kappa->Set(i, j, -(dnx_dx + dny_dy));
        }
    }

    // Face U: interface_u(i,j) btw phi(i-1,j) & phi(i,j)
    // phi(i-1, j) <-> Label(i,   j+1)
    // phi(i,   j) <-> Label(i+1, j+1)
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {   

            const labeltype left  = fields->Label(i,     j + 1);
            const labeltype right = fields->Label(i + 1, j + 1);

            const bool leftFluid  = IS_FLUID(left);
            const bool rightFluid = IS_FLUID(right);

            if (leftFluid == rightFluid) continue;
            if (IS_SOLID(left)) continue;
            if (IS_SOLID(right)) continue;

            const varType kappa_face = leftFluid ? fields->kappa->Get(i - 1, j)  
                                                 : fields->kappa->Get(i, j); 
            fields->interface_u->Set(i, j, -gamma * kappa_face);
        }
    }

    // faces V: interface_v(i,j) btw phi(i,j-1) & phi(i,j)
    // phi(i, j-1) <-> Label(i+1, j)
    // phi(i, j) <-> Label(i+1, j+1)
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 1; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const labeltype labelBot = fields->Label(i + 1, j);
            const labeltype labelTop = fields->Label(i + 1, j + 1);

            const bool botFluid = IS_FLUID(labelBot);
            const bool topFluid = IS_FLUID(labelTop);

            if (botFluid == topFluid) continue;
            if (IS_SOLID(labelBot))   continue;
            if (IS_SOLID(labelTop))   continue;

            const varType kappa_face = botFluid ? fields->kappa->Get(i, j - 1)   
                                                : fields->kappa->Get(i, j);      
            fields->interface_v->Set(i, j, -gamma * kappa_face);
        }
    }
}

void PIC::LaplacePressure() {
    UpdatePhiFromParticles();
    ComputeSurfaceTensionOnFaces();
}