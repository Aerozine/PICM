#include "SemiLagrangian.hpp"
#include <algorithm>
#include <cmath>

// Semi-Lagrangian advection
//  Each velocity component is advected independently:
//    1. For every face (i,j), trace a particle backward in time using RK2
//       to find the "departure point" (x_dep, y_dep).
//    2. Interpolate the current velocity field at that point.
//    3. Store the result in new grids, then move them into the fields.
//
//  Using separate new grids ensures all reads come from the current-step
//  values — equivalent to a Jacobi-style update.
//
//  Loop order: j (outer) → i (inner) so that consecutive Set() calls write
//  to consecutive memory locations (row-major: A[nx*j + i]).

void SemiLagrangian::Advect() const {
  Grid2D uNew(fields->u.nx, fields->u.ny);
  Grid2D vNew(fields->v.nx, fields->v.ny);

  for (int j = 0; j < fields->u.ny; ++j)
    for (int i = 0; i < fields->u.nx; ++i) {
      varType x, y;
      traceParticleU(i, j, x, y);
      uNew.Set(i, j, interpolateU(x, y));
    }

  for (int j = 0; j < fields->v.ny; ++j)
    for (int i = 0; i < fields->v.nx; ++i) {
      varType x, y;
      traceParticleV(i, j, x, y);
      vNew.Set(i, j, interpolateV(x, y));
    }

  fields->u = std::move(uNew);
  fields->v = std::move(vNew);
}

void SemiLagrangian::AdvectSmoke() const {
  Grid2D smokeNew(fields->smokeMap.nx, fields->smokeMap.ny);

  for (int j = 0; j < fields->smokeMap.ny; ++j) {
    for (int i = 0; i < fields->smokeMap.nx; ++i) {

      // Physical position of cell centre (i, j)
      const varType x0 = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
      const varType y0 = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;

      // RK2 backward trace
      varType u0, v0;
      getVelocity(x0, y0, u0, v0);
      const varType xMid = x0 - REAL_LITERAL(0.5) * dt * u0;
      const varType yMid = y0 - REAL_LITERAL(0.5) * dt * v0;

      varType uMid, vMid;
      getVelocity(xMid, yMid, uMid, vMid);
      varType xDep = x0 - dt * uMid;
      varType yDep = y0 - dt * vMid;

      xDep = std::clamp(xDep, REAL_LITERAL(0.0),
                        static_cast<varType>(nx - 1) * dx);
      yDep = std::clamp(yDep, REAL_LITERAL(0.0),
                        static_cast<varType>(ny - 1) * dy);

      smokeNew.Set(i, j, interpolateSmoke(xDep, yDep));
    }
  }

  fields->smokeMap = std::move(smokeNew);
}

// RK2 backward particle traces

void SemiLagrangian::traceParticleU(const int i, const int j, varType &x,
                                    varType &y) const {
  // u-face physical position: (i·dx, (j+0.5)·dy).
  const varType x0 = static_cast<varType>(i) * dx;
  const varType y0 = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;

  varType u0, v0;
  getVelocity(x0, y0, u0, v0);
  const varType xMid = x0 - REAL_LITERAL(0.5) * dt * u0;
  const varType yMid = y0 - REAL_LITERAL(0.5) * dt * v0;

  varType uMid, vMid;
  getVelocity(xMid, yMid, uMid, vMid);
  x = x0 - dt * uMid;
  y = y0 - dt * vMid;

  x = std::clamp(x, REAL_LITERAL(0.0), static_cast<varType>(nx - 1) * dx);
  y = std::clamp(y, REAL_LITERAL(0.0), static_cast<varType>(ny - 1) * dy);
}

void SemiLagrangian::traceParticleV(const int i, const int j, varType &x,
                                    varType &y) const {
  // v-face physical position: ((i+0.5)·dx, j·dy).
  const varType x0 = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
  const varType y0 = static_cast<varType>(j) * dy;

  varType u0, v0;
  getVelocity(x0, y0, u0, v0);
  const varType xMid = x0 - REAL_LITERAL(0.5) * dt * u0;
  const varType yMid = y0 - REAL_LITERAL(0.5) * dt * v0;

  varType uMid, vMid;
  getVelocity(xMid, yMid, uMid, vMid);
  x = x0 - dt * uMid;
  y = y0 - dt * vMid;

  x = std::clamp(x, REAL_LITERAL(0.0), static_cast<varType>(nx - 1) * dx);
  y = std::clamp(y, REAL_LITERAL(0.0), static_cast<varType>(ny - 1) * dy);
}

// Bilinear interpolation

varType SemiLagrangian::interpolateU(const varType x, const varType y) const {
  const varType i_real = x / dx;
  const varType j_real = y / dy - REAL_LITERAL(0.5);

  int i = static_cast<int>(std::floor(i_real));
  int j = static_cast<int>(std::floor(j_real));

  const varType fx = i_real - static_cast<varType>(i);
  const varType fy = j_real - static_cast<varType>(j);

  i = std::clamp(i, 0, fields->u.nx - 2);
  j = std::clamp(j, 0, fields->u.ny - 2);

  const varType u00 = fields->u.Get(i, j);
  const varType u10 = fields->u.Get(i + 1, j);
  const varType u01 = fields->u.Get(i, j + 1);
  const varType u11 = fields->u.Get(i + 1, j + 1);

  return (REAL_LITERAL(1.0) - fy) *
             ((REAL_LITERAL(1.0) - fx) * u00 + fx * u10) +
         fy * ((REAL_LITERAL(1.0) - fx) * u01 + fx * u11);
}

varType SemiLagrangian::interpolateV(const varType x, const varType y) const {
  const varType i_real = x / dx - REAL_LITERAL(0.5);
  const varType j_real = y / dy;

  int i = static_cast<int>(std::floor(i_real));
  int j = static_cast<int>(std::floor(j_real));

  const varType fx = i_real - static_cast<varType>(i);
  const varType fy = j_real - static_cast<varType>(j);

  i = std::clamp(i, 0, fields->v.nx - 2);
  j = std::clamp(j, 0, fields->v.ny - 2);

  const varType v00 = fields->v.Get(i, j);
  const varType v10 = fields->v.Get(i + 1, j);
  const varType v01 = fields->v.Get(i, j + 1);
  const varType v11 = fields->v.Get(i + 1, j + 1);

  return (REAL_LITERAL(1.0) - fy) *
             ((REAL_LITERAL(1.0) - fx) * v00 + fx * v10) +
         fy * ((REAL_LITERAL(1.0) - fx) * v01 + fx * v11);
}

void SemiLagrangian::getVelocity(const varType x, const varType y, varType &u,
                                 varType &v) const {
  u = interpolateU(x, y);
  v = interpolateV(x, y);
}

varType SemiLagrangian::interpolateSmoke(const varType x,
                                         const varType y) const {
  // smokeMap is cell-centred: (i+0.5)*dx, (j+0.5)*dy
  const varType i_real = x / dx - REAL_LITERAL(0.5);
  const varType j_real = y / dy - REAL_LITERAL(0.5);

  int i = static_cast<int>(std::floor(i_real));
  int j = static_cast<int>(std::floor(j_real));

  const varType fx = i_real - static_cast<varType>(i);
  const varType fy = j_real - static_cast<varType>(j);

  i = std::clamp(i, 0, fields->smokeMap.nx - 2);
  j = std::clamp(j, 0, fields->smokeMap.ny - 2);

  const varType s00 = fields->smokeMap.Get(i, j);
  const varType s10 = fields->smokeMap.Get(i + 1, j);
  const varType s01 = fields->smokeMap.Get(i, j + 1);
  const varType s11 = fields->smokeMap.Get(i + 1, j + 1);

  return (REAL_LITERAL(1.0) - fy) *
             ((REAL_LITERAL(1.0) - fx) * s00 + fx * s10) +
         fy * ((REAL_LITERAL(1.0) - fx) * s01 + fx * s11);
}
