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


    fields->phi->reset();
    fields->kappa->reset();
    fields->normalX->reset();
    fields->normalY->reset();

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

// Paramètre à ajouter dans params (ou en constante pour tester)
// params.contactAngle en radians  (ex: M_PI/4 = 45°, M_PI/2 = 90°, 2*M_PI/3 = 120°)

void PIC::ComputeSurfaceTensionOnFaces() const
{
    const int nx = fields->nx;
    const int ny = fields->ny;
    const varType dx = fields->dx;
    const varType dy = fields->dy;
    const varType gamma = params.gamma;
    const varType cosTheta = std::cos(params.contactAngle); // NOUVEAU

    fields->interface_u->reset();
    fields->interface_v->reset();

    // ----------------------------------------------------------------
    // Helpers : phi fantôme au bord solide
    //
    // Si le voisin est solide, on ne clampe plus phi(voisin) = phi(i,j)
    // On utilise à la place :  phi_ghost = phi(i,j) - d * cosTheta
    // où d est le pas dans la direction normale au solide (dx ou dy).
    //
    // Cela impose ∂φ/∂n_solide = -cosTheta au premier ordre.
    // ----------------------------------------------------------------

    // phi fantôme en X : voisin (i+di, j), solide ou hors domaine
    auto phiGhostX = [&](int i, int j, int di) -> varType {
        int ni = i + di;
        // hors domaine
        if (ni < 0 || ni >= nx)
            return fields->phi->Get(i, j) - di * dx * cosTheta;
        // voisin solide
        if (IS_SOLID(fields->Label(ni + 1, j + 1)))
            return fields->phi->Get(i, j) - di * dx * cosTheta;
        // cas normal
        return fields->phi->Get(ni, j);
    };

    // phi fantôme en Y : voisin (i, j+dj), solide ou hors domaine
    auto phiGhostY = [&](int i, int j, int dj) -> varType {
        int nj = j + dj;
        if (nj < 0 || nj >= ny)
            return fields->phi->Get(i, j) - dj * dy * cosTheta;
        if (IS_SOLID(fields->Label(i + 1, nj + 1)))
            return fields->phi->Get(i, j) - dj * dy * cosTheta;
        return fields->phi->Get(i, nj);
    };

    // ----------------------------------------------------------------
    // Normales — même structure qu'avant, mais via les helpers
    // ----------------------------------------------------------------
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {

            const varType phi_xp = phiGhostX(i, j, +1); // MODIFIÉ
            const varType phi_xm = phiGhostX(i, j, -1); // MODIFIÉ
            const varType phi_yp = phiGhostY(i, j, +1); // MODIFIÉ
            const varType phi_ym = phiGhostY(i, j, -1); // MODIFIÉ

            const varType gx = (phi_xp - phi_xm) / (2.0 * dx);
            const varType gy = (phi_yp - phi_ym) / (2.0 * dy);
            const varType gn = std::hypot(gx, gy);

            if (gn > REAL_EPSILON) {
                fields->normalX->Set(i, j, gx / gn);
                fields->normalY->Set(i, j, gy / gn);
            }
        }
    }

    // ----------------------------------------------------------------
    // Courbure — idem, on skip les cellules solides comme avant
    // ----------------------------------------------------------------
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {

            if (!IS_FLUID(fields->Label(i + 1, j + 1))) continue;

            // Pour la divergence des normales on utilise aussi les fantômes,
            // mais on travaille sur normalX/normalY déjà calculés.
            // Pour les voisins solides de normalX/normalY, on n'a pas de valeur :
            // on clampe à 0 (normal tangente au solide → divergence nulle là).
            auto nxVal = [&](int ii, int jj) -> varType {
                if (ii < 0 || ii >= nx) return 0.0;
                if (IS_SOLID(fields->Label(ii + 1, jj + 1))) return 0.0;
                return fields->normalX->Get(ii, jj);
            };
            auto nyVal = [&](int ii, int jj) -> varType {
                if (jj < 0 || jj >= ny) return 0.0;
                if (IS_SOLID(fields->Label(ii + 1, jj + 1))) return 0.0;
                return fields->normalY->Get(ii, jj);
            };

            const varType dnx_dx =
                (nxVal(i + 1, j) - nxVal(i - 1, j)) / (2.0 * dx);
            const varType dny_dy =
                (nyVal(i, j + 1) - nyVal(i, j - 1)) / (2.0 * dy);

            fields->kappa->Set(i, j, -(dnx_dx + dny_dy));
        }
    }

    // ----------------------------------------------------------------
    // Faces U et V — inchangées
    // ----------------------------------------------------------------
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            const labeltype left  = fields->Label(i,     j + 1);
            const labeltype right = fields->Label(i + 1, j + 1);
            const bool leftFluid  = IS_FLUID(left);
            const bool rightFluid = IS_FLUID(right);
            if (leftFluid == rightFluid) continue;
            if (IS_SOLID(left))  continue;
            if (IS_SOLID(right)) continue;
            const varType kappa_face = leftFluid ? fields->kappa->Get(i - 1, j)
                                                 : fields->kappa->Get(i,     j);
            fields->interface_u->Set(i, j, -gamma * kappa_face);
        }
    }

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