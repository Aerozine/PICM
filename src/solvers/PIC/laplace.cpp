#include "PIC.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
/*
// ============================================================
//  Surface tension CSF — Brackbill et al. (1992)
//
//  Δp = γ · κ      (Young-Laplace)
//
//  Méthode :
//    1. Calculer une color function C(i,j) = densité de
//       particules lissée par un kernel SPH sur la grille
//    2. Normaliser C ∈ [0,1]  (0=air, 1=fluide plein)
//    3. Calculer le gradient ∇C par différences finies centrées
//    4. Normaliser : n = ∇C / |∇C|  (normale à l'interface)
//    5. Courbure : κ = -∇·n
//    6. Force volumique : f = γ · κ · ∇C
//    7. Distribuer f aux particules via interpolation bilinéaire
// ============================================================

static const varType CSF_KERNEL_RADIUS_FACTOR = 2.0; // h = factor * dx

// Kernel de lissage gaussien tronqué pour la color function
static inline varType smoothingKernel(varType r, varType h)
{
    const varType q = r / h;
    if (q >= 1.0) return 0.0;
    // Kernel cubique (Monaghan) normalisé 2D
    const varType s = 1.0 - q;
    return s * s * s;   // normalisation absorbée dans la somme
}

void PIC::particleInteraction()
{
    const int    nx = fields->nx;
    const int    ny = fields->ny;
    const varType dx = fields->dx;
    const varType dy = fields->dy;   // supposé = dx pour grille carrée

    const varType gamma = 0.073;  // [N/m]
    const varType mp    = 0.0000001;
    const varType h_smooth = CSF_KERNEL_RADIUS_FACTOR * dx;
    const int     Rsmooth  = static_cast<int>(std::ceil(h_smooth / dx)) + 1;

    const varType maxDeltaV = 0.2 * dx;

    // --------------------------------------------------------
    // Étape 1 : construire la color function C sur la grille
    //   C(i,j) = somme des kernels de chaque particule
    //            évaluée au centre de la cellule (i,j)
    // --------------------------------------------------------
    // Indexation linéaire : idx(i,j) = (j-1)*nx + (i-1)
    const int N = (nx + 2) * (ny + 2);   // avec ghost cells
    auto IDX = [&](int i, int j) { return j * (nx + 2) + i; };

    std::vector<varType> C(N, 0.0);
    std::vector<varType> Cx(N, 0.0);  // ∂C/∂x
    std::vector<varType> Cy(N, 0.0);  // ∂C/∂y
    std::vector<varType> kappa(N, 0.0);

    // Accumulation de la color function
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {

            if (!IS_FLUID(fields->Label(i, j))) continue;

            const varType xc = (i - 0.5) * dx;  // centre cellule
            const varType yc = (j - 0.5) * dy;

            const int ni_lo = std::max(1, i - Rsmooth);
            const int ni_hi = std::min(nx, i + Rsmooth);
            const int nj_lo = std::max(1, j - Rsmooth);
            const int nj_hi = std::min(ny, j + Rsmooth);

            varType cval = 0.0;

            for (int nj = nj_lo; nj <= nj_hi; ++nj) {
                for (int ni = ni_lo; ni <= ni_hi; ++ni) {

                    const Particles &cell = (*cloud)(ni, nj);
                    const int np = cell.size();

                    for (int p = 0; p < np; ++p) {
                        const varType xp = cell.GetX(p);
                        const varType yp = cell.GetY(p);
                        const varType r  = std::hypot(xp - xc, yp - yc);
                        cval += smoothingKernel(r, h_smooth);
                    }
                }
            }

            C[IDX(i, j)] = cval;
        }
    }

    // Normaliser C dans [0, 1]
    varType Cmax = *std::max_element(C.begin(), C.end());
    if (Cmax < REAL_EPSILON) return;   // pas de fluide
    OMP_PRAGMA(omp parallel for schedule(static))
    for (int k = 0; k < N; ++k)
        C[k] /= Cmax;

    // --------------------------------------------------------
    // Étape 2 : gradient de C par différences finies centrées
    //   ∂C/∂x ≈ (C(i+1,j) - C(i-1,j)) / (2·dx)
    // --------------------------------------------------------
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            Cx[IDX(i,j)] = (C[IDX(i+1,j)] - C[IDX(i-1,j)]) / (2.0 * dx);
            Cy[IDX(i,j)] = (C[IDX(i,j+1)] - C[IDX(i,j-1)]) / (2.0 * dy);
        }
    }

    // --------------------------------------------------------
    // Étape 3 : courbure κ = -∇·(∇C / |∇C|)
    //
    //   n = ∇C / |∇C|
    //   κ = -(∂nx/∂x + ∂ny/∂y)
    //
    //   On calcule d'abord nx,ny partout,
    //   puis on différencie.
    // --------------------------------------------------------
    std::vector<varType> nx_field(N, 0.0);
    std::vector<varType> ny_field(N, 0.0);

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            const varType gx = Cx[IDX(i,j)];
            const varType gy = Cy[IDX(i,j)];
            const varType gn = std::hypot(gx, gy);
            if (gn > REAL_EPSILON) {
                nx_field[IDX(i,j)] = gx / gn;
                ny_field[IDX(i,j)] = gy / gn;
            }
        }
    }

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            if (!IS_FLUID(fields->Label(i,j))) continue;

            const varType dnx_dx = (nx_field[IDX(i+1,j)] - nx_field[IDX(i-1,j)]) / (2.0 * dx);
            const varType dny_dy = (ny_field[IDX(i,j+1)] - ny_field[IDX(i,j-1)]) / (2.0 * dy);

            kappa[IDX(i,j)] = -(dnx_dx + dny_dy);
        }
    }

    // --------------------------------------------------------
    // Étape 4 : force CSF sur les particules
    //
    //   f_vol = γ · κ · ∇C      [N/m³]
    //
    //   Pour convertir en Δv sur une particule :
    //   Δv = f_vol / ρ_local * dt
    //
    //   On approche ρ_local par le nombre de particules
    //   dans la cellule × mp / (dx·dy).
    //   On distribue f_vol à chaque particule par
    //   interpolation bilinéaire depuis la grille.
    // --------------------------------------------------------
    OMP_PRAGMA(omp parallel for collapse(2) schedule(dynamic, 1))
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {

            if (!IS_FLUID(fields->Label(i,j))) continue;

            Particles &cell = (*cloud)(i - 1, j -1);
            const int np = cell.size();
            if (np == 0) continue;

            // Densité locale estimée
            const varType rho_local = np * mp / (dx * dy);
            if (rho_local < REAL_EPSILON) continue;

            for (int p = 0; p < np; ++p) {

                const varType xp = cell.GetX(p);
                const varType yp = cell.GetY(p);

                // Position normalisée dans la cellule pour interp bilinéaire
                // On interpole kappa et ∇C aux 4 coins de la cellule
                // Ici simplification : on utilise la valeur au centre (ordre 1)
                // Pour ordre 2 : interpoler depuis les 4 cellules voisines

                // Interpolation bilinéaire de κ et ∇C à la position (xp, yp)
                // Indices des 4 cellules voisines encadrantes
                const int il = std::max(1, static_cast<int>(std::floor(xp / dx)));
                const int ir = std::min(nx, il + 1);
                const int jd = std::max(1, static_cast<int>(std::floor(yp / dy)));
                const int ju = std::min(ny, jd + 1);

                const varType tx = (xp - (il - 0.5) * dx) / dx;
                const varType ty = (yp - (jd - 0.5) * dy) / dy;
                const varType tx1 = 1.0 - tx, ty1 = 1.0 - ty;

                // Poids bilinéaires
                const varType w00 = tx1 * ty1;
                const varType w10 = tx  * ty1;
                const varType w01 = tx1 * ty;
                const varType w11 = tx  * ty;

                // Interpoler κ
                const varType kap =
                    w00 * kappa[IDX(il, jd)] + w10 * kappa[IDX(ir, jd)] +
                    w01 * kappa[IDX(il, ju)] + w11 * kappa[IDX(ir, ju)];

                // Interpoler ∇C
                const varType gcx =
                    w00 * Cx[IDX(il, jd)] + w10 * Cx[IDX(ir, jd)] +
                    w01 * Cx[IDX(il, ju)] + w11 * Cx[IDX(ir, ju)];
                const varType gcy =
                    w00 * Cy[IDX(il, jd)] + w10 * Cy[IDX(ir, jd)] +
                    w01 * Cy[IDX(il, ju)] + w11 * Cy[IDX(ir, ju)];

                // Force volumique CSF : f = γ · κ · ∇C
                const varType fx = gamma * kap * gcx;
                const varType fy = gamma * kap * gcy;

                // Δv = f / ρ  (par unité de temps, supposé dt=1 intégré ailleurs)
                const varType du = std::clamp(fx / rho_local, -maxDeltaV, maxDeltaV);
                const varType dv = std::clamp(fy / rho_local, -maxDeltaV, maxDeltaV);

                cell.SetU(p, cell.GetU(p) + du);
                cell.SetV(p, cell.GetV(p) + dv);
            }
        }
    }
}
*/
