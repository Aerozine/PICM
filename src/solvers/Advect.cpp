#include "Solver.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>

// TODO: advect in SOLIDS is useless | add if(SOLID) {skip} ?
// is branching worse than looking in each solid ?
void Solver::Advect() {
  Grid2D uNew(fields->u.nx, fields->u.ny);
  Grid2D vNew(fields->v.nx, fields->v.ny);

  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < fields->u.ny; ++j)
    for (int i = 0; i < fields->u.nx; ++i) {
      // if label (i + 1, j + 1) is a BC_U, this means there is a
      // BC on the left boundary of this cell i.e. u(i, j)
      if (IS_BC_U(fields->Label(i, j + 1))) {
        uNew.Set(i, j, fields->u.Get(i, j));
        continue;
      }
      if (IS_SOLID(fields->Label(i, j + 1)) ||
          IS_SOLID(fields->Label(i + 1, j + 1))) {
        uNew.Set(i, j, 0.0);
        continue;
      }
      varType x, y;
      traceParticle(i, j, x, y,0);
      uNew.Set(i, j,
      fields->u.interpolate<0>(x, y, dx, dy));
    }

  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < fields->v.ny; ++j)
    for (int i = 0; i < fields->v.nx; ++i) {
      if (IS_BC_V(fields->Label(i + 1, j))) {
        assert(std::isfinite(fields->v.Get(i, j)));
        vNew.Set(i, j, fields->v.Get(i, j));
        continue;
      }
      if (IS_SOLID(fields->Label(i + 1, j)) ||
          IS_SOLID(fields->Label(i + 1, j + 1))) {
        vNew.Set(i, j, 0.0);
        continue;
      }

      varType x, y;
      traceParticle(i, j, x, y,1);
      assert(std::isfinite(fields->v.interpolate<1>(x, y, dx, dy)));
      vNew.Set(i, j,
      fields->v.interpolate<1>(x, y, dx, dy));
    }

  fields->u = std::move(uNew);
  fields->v = std::move(vNew);
}



// RK2 backward particle traces
// field: 0 for u-face, 1 for v-face
void Solver::traceParticle(const int i, const int j, varType &x, varType &y,
                           const uint8_t field) const {
  varType x0, y0;

  if (field == 0) {
    // u-face physical position: (i·dx, (j+0.5)·dy)
    x0 = static_cast<varType>(i) * dx;
    y0 = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;
  } else {
    // v-face physical position: ((i+0.5)·dx, j·dy)
    x0 = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
    y0 = static_cast<varType>(j) * dy;
  }

  // RK2 backward trace: first half-step
  varType u0, v0;
  getVelocity(x0, y0, u0, v0);
  const varType xMid = x0 - REAL_LITERAL(0.5) * dt * u0;
  const varType yMid = y0 - REAL_LITERAL(0.5) * dt * v0;

  // RK2 backward trace: full step with midpoint velocity
  varType uMid, vMid;
  getVelocity(xMid, yMid, uMid, vMid);
  x = x0 - dt * uMid;
  y = y0 - dt * vMid;

  // Clamp to valid range
  // @todo for a particle u should not be nx-1,ny-1 the size is not the same
  // 1 cell is thrown away
  if (field == 0) {
    x = std::clamp(x, REAL_LITERAL(0.0), static_cast<varType>(fields->u.nx) * dx);
    y = std::clamp(y, REAL_LITERAL(0.0), static_cast<varType>(fields->u.ny) * dy);
  } else {
    x = std::clamp(x, REAL_LITERAL(0.0), static_cast<varType>(fields->v.nx) * dx);
    y = std::clamp(y, REAL_LITERAL(0.0), static_cast<varType>(fields->v.ny) * dy);
  }
}
