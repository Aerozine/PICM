#include "PIC.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

// we can look at the cloud to see if there is particle
// it there is we can simply set fluid
void PIC::UpdateCellState() const {
  // In a fully-filled simulation the initial labels are the domain state.
  // Do not turn temporarily empty PIC cells into AIR: refill may add
  // particles after this pass, and the pressure solve would then see a
  // one-step stale air column at inflows.
  if (!params.freeSurface)
    return;

    OMP_PRAGMA(omp parallel for collapse(2))
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        // early exit
        const labeltype current = fields->Label(i + 1, j + 1);
        if (IS_SOLID(current) || IS_BC_U(current) || IS_BC_V(current))
          continue;
        if (cloud->countIn(i, j) > 0) {
          fields->setFluid(i + 1, j + 1);
          continue;
        }
        fields->setAir(i + 1, j + 1);
      }
    }

    // trick to help contact angles
    if (params.surfaceTension) {
      const int targetPPC = std::max(1, params.ppcx * params.ppcy);
      const int sparseWallThreshold = std::max(1, targetPPC / 4);
      std::vector<unsigned char> removeWallSpeck(
          static_cast<std::size_t>(nx) * ny, 0);

      const auto cellIndex = [this](int i, int j) {
        return static_cast<std::size_t>(nx) * j + i;
      };
      const auto isSolidCell = [&](int i, int j) {
        if (i < 0 || i >= nx || j < 0 || j >= ny)
          return true;
        return IS_SOLID(fields->Label(i + 1, j + 1));
      };
      const auto isFluidCell = [&](int i, int j) {
        if (i < 0 || i >= nx || j < 0 || j >= ny)
          return false;
        return IS_FLUID(fields->Label(i + 1, j + 1));
      };

      OMP_PRAGMA(omp parallel for collapse(2))
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          if (!isFluidCell(i, j) || cloud->countIn(i, j) >= sparseWallThreshold)
            continue;

          bool touchesSolid = false;
          bool supportedAwayFromWall = false;
          const int dirs[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
          for (const auto &dir : dirs) {
            const int si = i + dir[0];
            const int sj = j + dir[1];
            if (!isSolidCell(si, sj))
              continue;

            touchesSolid = true;
            const int oi = i - dir[0];
            const int oj = j - dir[1];
            if (isFluidCell(oi, oj) &&
                cloud->countIn(oi, oj) >= sparseWallThreshold) {
              supportedAwayFromWall = true;
              break;
            }
          }

          if (touchesSolid && !supportedAwayFromWall)
            removeWallSpeck[cellIndex(i, j)] = 1;
        }
      }

      OMP_PRAGMA(omp parallel for collapse(2))
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
          if (removeWallSpeck[cellIndex(i, j)])
            fields->setAir(i + 1, j + 1);
    }

    OMP_PRAGMA(omp parallel for collapse(2))
    for (int i = 1; i < fields->p.nx - 1; ++i) {
      for (int j = 1; j < fields->p.ny - 1; ++j) {
        if (!IS_AIR(fields->Label(i, j)))
          continue;

        fields->p.Set(i, j, 0.0);

        labeltype leftLabel = fields->Label(i - 1, j);
        labeltype rightLabel = fields->Label(i + 1, j);
        labeltype bottomLabel = fields->Label(i, j - 1);
        labeltype topLabel = fields->Label(i, j + 1);

        // u, v = 0 at AIR|AIR or AIR|SOLID interfaces
        if (!IS_FLUID(leftLabel))
          fields->u.Set(i - 1, j - 1, static_cast<varType>(0));
        if (!IS_FLUID(rightLabel))
          fields->u.Set(i, j - 1, static_cast<varType>(0));
        if (!IS_FLUID(bottomLabel))
          fields->v.Set(i - 1, j - 1, static_cast<varType>(0));
        if (!IS_FLUID(topLabel))
          fields->v.Set(i - 1, j, static_cast<varType>(0));
      }
    }
}
