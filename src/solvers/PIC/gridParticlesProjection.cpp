#include "PIC.hpp"
#include <iostream>

// ScatterToGrid kept for external callers (FLIP, APIC subclasses).
// No longer called by ProjectParticlesOnGrid in PIC — the 4-color cell
// loop below eliminates the need for atomics entirely.
void PIC::ScatterToGrid(varType xg, varType yg, varType val, Grid2D &sum,
                        Grid2D &weight, int imax, int jmax) {
    const int i0 = static_cast<int>(std::floor(xg));
    const int j0 = static_cast<int>(std::floor(yg));

#ifdef USE_SPEED
    constexpr int8_t radius = 1;
#else
    constexpr int8_t radius = 2;
#endif

    for (int8_t dj = -radius; dj <= radius; ++dj) {
        for (int8_t di = -radius; di <= radius; ++di) {
            const int i = i0 + di, j = j0 + dj;
            if (i < 0 || i >= imax || j < 0 || j >= jmax)
                continue;
            const varType k = hat(xg - varType(i)) * hat(yg - varType(j));
            if (k > varType(0)) {
                OMP_PRAGMA(omp atomic)
                sum(i, j) += k * val;
                OMP_PRAGMA(omp atomic)
                weight(i, j) += k;
            }
        }
    }
}

// Scatter all particles in one cell into the MAC grid.
// Uses atomics because the kernel radius (1 or 2) means that particles in
// neighbouring cells write to overlapping grid nodes.  Specifically, a
// particle in cell ci reaches u-nodes [ci-R, ci+R] and the v y-stagger
// shifts the j-footprint by an extra ±1, so a safe conflict-free coloring
// would need stride ≥ 2R+2 in both dimensions (≥6 for R=2, ≥4 for R=1).
// A plain parallel-for with atomics is simpler and correct.
void PIC::scatterCell(const Particles &cell,
                      Grid2D &u_sum, Grid2D &u_weight,
                      Grid2D &v_sum, Grid2D &v_weight) {
#ifdef USE_SPEED
    constexpr int8_t radius = 1;
#else
    constexpr int8_t radius = 2;
#endif

    for (int p = 0; p < cell.size(); ++p) {
        const varType x  = cell.GetX(p);
        const varType y  = cell.GetY(p);
        const varType up = cell.GetU(p);
        const varType vp = cell.GetV(p);

        // u-face stagger: xg = x/dx,  yg = y/dy - 0.5
        {
            const varType xg = x / dx;
            const varType yg = y / dy - varType(0.5);
            const int i0 = static_cast<int>(std::floor(xg));
            const int j0 = static_cast<int>(std::floor(yg));
            for (int8_t dj = -radius; dj <= radius; ++dj) {
                for (int8_t di = -radius; di <= radius; ++di) {
                    const int i = i0 + di, j = j0 + dj;
                    if (i < 0 || i >= nx + 1 || j < 0 || j >= ny) continue;
                    const varType k = hat(xg - varType(i)) * hat(yg - varType(j));
                    if (k > varType(0)) {
                        OMP_PRAGMA(omp atomic)
                        u_sum(i, j)    += k * up;
                        OMP_PRAGMA(omp atomic)
                        u_weight(i, j) += k;
                    }
                }
            }
        }

        // v-face stagger: xg = x/dx - 0.5,  yg = y/dy
        {
            const varType xg = x / dx - varType(0.5);
            const varType yg = y / dy;
            const int i0 = static_cast<int>(std::floor(xg));
            const int j0 = static_cast<int>(std::floor(yg));
            for (int8_t dj = -radius; dj <= radius; ++dj) {
                for (int8_t di = -radius; di <= radius; ++di) {
                    const int i = i0 + di, j = j0 + dj;
                    if (i < 0 || i >= nx || j < 0 || j >= ny + 1) continue;
                    const varType k = hat(xg - varType(i)) * hat(yg - varType(j));
                    if (k > varType(0)) {
                        OMP_PRAGMA(omp atomic)
                        v_sum(i, j)    += k * vp;
                        OMP_PRAGMA(omp atomic)
                        v_weight(i, j) += k;
                    }
                }
            }
        }
    }
}

void PIC::ProjectParticlesOnGrid() {
    fields->u_sum->reset();
    fields->u_weight->reset();
    fields->v_sum->reset();
    fields->v_weight->reset();

    // Parallel scatter: each thread owns one cell and calls scatterCell,
    // which uses atomics on the shared grid arrays.  The old 4-color
    // scheme attempted to avoid atomics, but the required conflict-free
    // stride is ≥2R+2 per dimension (≥6 for R=2, ≥4 for R=1), making
    // it far more than 4 colors — so atomics are the correct approach here.
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int ci = 0; ci < nx; ++ci) {
        for (int cj = 0; cj < ny; ++cj) {
            scatterCell((*cloud)(ci, cj),
                        *fields->u_sum, *fields->u_weight,
                        *fields->v_sum, *fields->v_weight);
        }
    }

    // normalise u
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < fields->u.ny; ++j) {
        for (int i = 0; i < fields->u.nx; ++i) {
            const labeltype left  = fields->Label(i,     j + 1);
            const labeltype right = fields->Label(i + 1, j + 1);

            if (IS_SOLID(left) || IS_SOLID(right)) {
                fields->u.Set(i, j, varType(0));
                continue;
            }
            if (IS_BC_U(left)) continue;

            const varType w = fields->u_weight->Get(i, j);
            fields->u.Set(i, j, w >= REAL_EPSILON
                                ? fields->u_sum->Get(i, j) / w
                                : varType(0));
        }
    }

    // normalise v
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < fields->v.ny; ++j) {
        for (int i = 0; i < fields->v.nx; ++i) {
            const labeltype bottom = fields->Label(i + 1, j);
            const labeltype top    = fields->Label(i + 1, j + 1);

            if (IS_SOLID(bottom) || IS_SOLID(top)) {
                fields->v.Set(i, j, varType(0));
                continue;
            }
            if (IS_BC_V(bottom)) continue;

            const varType w = fields->v_weight->Get(i, j);
            fields->v.Set(i, j, w >= REAL_EPSILON
                                ? fields->v_sum->Get(i, j) / w
                                : varType(0));
        }
    }
}

void PIC::ProjectGridOnParticles() {
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int ci = 0; ci < nx; ++ci) {
        for (int cj = 0; cj < ny; ++cj) {
            Particles &cell = (*cloud)(ci, cj);
            for (int p = 0; p < cell.size(); ++p) {
                const varType x = cell.GetX(p);
                const varType y = cell.GetY(p);
                cell.SetU(p, fields->u.interpolate<0>(x, y, dx, dy));
                cell.SetV(p, fields->v.interpolate<1>(x, y, dx, dy) - dt * params.gravity);
            }
        }
    }
}
