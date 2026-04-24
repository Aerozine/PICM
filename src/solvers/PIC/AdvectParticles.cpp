#include "PIC.hpp"
#include <cmath>

void PIC::Advect() {
    struct PendingParticle {
        varType x;
        varType y;
        varType u;
        varType v;
    };

    const varType xMax = dx * nx;
    const varType yMax = dy * ny;
    std::vector<std::vector<PendingParticle>> incoming(
        static_cast<std::size_t>(nx) * ny);

    // Each thread owns its source cell exclusively.
    // Cross-cell moves are staged into per-destination buffers, then merged
    // after the parallel loop so we never mutate a destination cell while
    // another thread is iterating over it.
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int ci = 0; ci < nx; ++ci) {
        for (int cj = 0; cj < ny; ++cj) {
            Particles &cell = (*cloud)(ci, cj);
            int local_idx = 0;
            int n = cell.size();

            while (local_idx < n) {
                const varType x0 = cell.GetX(local_idx);
                const varType y0 = cell.GetY(local_idx);
                const varType u0 = cell.GetU(local_idx);
                const varType v0 = cell.GetV(local_idx);

                // RK2 midpoint advection (same as before)
                const varType xmid = x0 + varType(0.5) * dt * u0;
                const varType ymid = y0 + varType(0.5) * dt * v0;

                const varType umid = fields->u.interpolate<0>(xmid, ymid, dx, dy);
                const varType vmid = fields->v.interpolate<1>(xmid, ymid, dx, dy);

                const varType x1 = x0 + dt * umid;
                const varType y1 = y0 + dt * vmid;

                // out of bounds: kill particle, swap-and-pop keeps array packed.
                // Decrement n — one fewer original particle remains.
                if (x1 < varType(0) || x1 >= xMax || y1 < varType(0) || y1 >= yMax) {
                    cell.Remove(local_idx);
                    --n;
                    continue;
                }

                const int ci1 = std::clamp(static_cast<int>(std::floor(x1 / dx)), 0, nx - 1);
                const int cj1 = std::clamp(static_cast<int>(std::floor(y1 / dy)), 0, ny - 1);

                if (IS_SOLID(fields->Label(ci1 + 1, cj1 + 1))) {
                    // stays at old position
                    ++local_idx;
                    continue;
                }

                if (ci1 != ci || cj1 != cj) {
                    const std::size_t dstIdx = static_cast<std::size_t>(ny) * ci1 + cj1;
#ifdef USE_OPENMP
                    omp_set_lock(&cloud->cellLocks[dstIdx]);
#endif
                    incoming[dstIdx].push_back(PendingParticle{x1, y1, u0, v0});
#ifdef USE_OPENMP
                    omp_unset_lock(&cloud->cellLocks[dstIdx]);
#endif

                    // particle left this source cell
                    cell.Remove(local_idx);
                    --n;
                    continue;
                }

                cell.SetX(local_idx, x1);
                cell.SetY(local_idx, y1);
                ++local_idx;
            }
        }
    }

    for (std::size_t dstIdx = 0; dstIdx < incoming.size(); ++dstIdx) {
        Particles &dst = cloud->cells[dstIdx];
        for (const PendingParticle &p : incoming[dstIdx])
            dst.Add(p.x, p.y, p.u, p.v, 0);
    }
}
