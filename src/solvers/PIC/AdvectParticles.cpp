#include "PIC.hpp"
#include <cmath>

void PIC::Advect() {
    const varType xMax = dx * nx;
    const varType yMax = dy * ny;

    // loop over cells — each thread owns its cell exclusively so
    // reads and in-place updates need no lock.
    // only cross-boundary transfers lock the destination cell.
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int ci = 0; ci < nx; ++ci) {
        for (int cj = 0; cj < ny; ++cj) {
            Particles &cell = (*cloud)(ci, cj);
            int local_idx = 0;
            // Capture the number of particles that belong to this cell at the
            // start of the timestep.  Other threads may call dst.Add() on this
            // cell while we iterate (transferring particles in), which raises
            // cell.size().  Without the cap those newly-arrived particles would
            // be re-advected a second time this step (double-advection bug).
            // We decrement n whenever a particle leaves the cell (Remove or
            // transfer) so the swap-and-pop bookkeeping stays correct.
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

                // update position in place before a possible transfer
                cell.SetX(local_idx, x1);
                cell.SetY(local_idx, y1);

                if (ci1 != ci || cj1 != cj) {
                    // particle crossed a cell boundary — transfer locks only dst.
                    // Decrement n — the particle has left this cell.
                    cloud->transfer(local_idx, ci, cj, ci1, cj1);
                    --n;
                    continue;
                }

                ++local_idx;
            }
        }
    }
    // countAliveParticles is now cloud->countIn(i,j) — no separate pass needed
}
