#include "PIC.hpp"
#include <cmath>

void PIC::Advect() {
  struct PendingParticle {
    varType x;
    varType y;
    varType u;
    varType v;
    varType cuX;
    varType cuY;
    varType cvX;
    varType cvY;
  };

  const varType xMax = dx * nx;
  const varType yMax = dy * ny;
  std::vector<std::vector<PendingParticle>> incoming(
      static_cast<std::size_t>(nx) * ny);
  const auto isFluidPoint = [&](const varType x, const varType y) {
    if (x < varType(0) || x >= xMax || y < varType(0) || y >= yMax)
      return false;
    const int ci = std::clamp(static_cast<int>(std::floor(x / dx)), 0, nx - 1);
    const int cj = std::clamp(static_cast<int>(std::floor(y / dy)), 0, ny - 1);
    return !IS_SOLID(fields->Label(ci + 1, cj + 1));
  };

  const auto clampToLastFluidPoint = [&](const varType x0, const varType y0,
                                         const varType x1, const varType y1,
                                         varType &xc, varType &yc) {
    varType lo = varType(0);
    varType hi = varType(1);
    for (int iter = 0; iter < 12; ++iter) {
      const varType mid = varType(0.5) * (lo + hi);
      const varType xm = x0 + mid * (x1 - x0);
      const varType ym = y0 + mid * (y1 - y0);
      if (isFluidPoint(xm, ym))
        lo = mid;
      else
        hi = mid;
    }

    const varType backedOff = std::max(varType(0), lo - REAL_LITERAL(1e-4));
    xc = x0 + backedOff * (x1 - x0);
    yc = y0 + backedOff * (y1 - y0);
    if (!isFluidPoint(xc, yc)) {
      xc = x0;
      yc = y0;
    }
  };

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
          const varType cuX0 = cell.GetCuX(local_idx);
          const varType cuY0 = cell.GetCuY(local_idx);
          const varType cvX0 = cell.GetCvX(local_idx);
          const varType cvY0 = cell.GetCvY(local_idx);

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

          varType xAdv = x1;
          varType yAdv = y1;
          if (!isFluidPoint(xAdv, yAdv))
            clampToLastFluidPoint(x0, y0, x1, y1, xAdv, yAdv);

          const int ci1 =
              std::clamp(static_cast<int>(std::floor(xAdv / dx)), 0, nx - 1);
          const int cj1 =
              std::clamp(static_cast<int>(std::floor(yAdv / dy)), 0, ny - 1);

          if (ci1 != ci || cj1 != cj) {
            const std::size_t dstIdx = static_cast<std::size_t>(ny) * ci1 + cj1;
#ifdef USE_OPENMP
            omp_set_lock(&cloud->cellLocks[dstIdx]);
#endif
            incoming[dstIdx].push_back(
                PendingParticle{xAdv, yAdv, u0, v0, cuX0, cuY0, cvX0, cvY0});
#ifdef USE_OPENMP
            omp_unset_lock(&cloud->cellLocks[dstIdx]);
#endif

            // particle left this source cell
            cell.Remove(local_idx);
            --n;
            continue;
          }

          cell.SetX(local_idx, xAdv);
          cell.SetY(local_idx, yAdv);
          ++local_idx;
        }
      }
    }

    for (std::size_t dstIdx = 0; dstIdx < incoming.size(); ++dstIdx) {
      Particles &dst = cloud->cells[dstIdx];
      for (const PendingParticle &p : incoming[dstIdx])
        dst.Add(p.x, p.y, p.u, p.v, 0, p.cuX, p.cuY, p.cvX, p.cvY);
    }
}
