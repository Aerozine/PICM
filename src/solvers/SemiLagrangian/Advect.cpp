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

void SemiLagrangian::Advect() const {
  Grid2D uNew(fields->u.nx, fields->u.ny);
  Grid2D vNew(fields->v.nx, fields->v.ny);

  for (int i = 0; i < fields->u.nx; ++i)
    for (int j = 0; j < fields->u.ny; ++j) {
      varType x, y;
      traceParticleU(i, j, x, y);
      uNew.Set(i, j, interpolateU(x, y));
    }

  for (int i = 0; i < fields->v.nx; ++i)
    for (int j = 0; j < fields->v.ny; ++j) {
      varType x, y;
      traceParticleV(i, j, x, y);
      vNew.Set(i, j, interpolateV(x, y));
    }

  fields->u = std::move(uNew);
  fields->v = std::move(vNew);
}

// RK2 backward particle traces

void SemiLagrangian::traceParticleU(const int i, const int j, varType &x,
                                    varType &y) const {
  // u-face physical position: (i.dx, (j+0.5).dy).
  const varType x0 = static_cast<varType>(i) * dx;
  const varType y0 = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;

  // Step 1: Euler half-step to midpoint.
  varType u0, v0;
  getVelocity(x0, y0, u0, v0);
  const varType xMid = x0 - REAL_LITERAL(0.5) * dt * u0;
  const varType yMid = y0 - REAL_LITERAL(0.5) * dt * v0;

  // Step 2: Full backward step using midpoint velocity.
  varType uMid, vMid;
  getVelocity(xMid, yMid, uMid, vMid);
  x = x0 - dt * uMid;
  y = y0 - dt * vMid;

  // Clamp to the physical domain so interpolation indices stay valid.
  x = std::clamp(x, REAL_LITERAL(0.0), static_cast<varType>(nx - 1) * dx);
  y = std::clamp(y, REAL_LITERAL(0.0), static_cast<varType>(ny - 1) * dy);
}

void SemiLagrangian::traceParticleV(const int i, const int j, varType &x,
                                    varType &y) const {
  // v-face physical position: ((i+0.5).dx, j.dy).
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
//
// Both functions follow the same pattern:
//   1. Map physical (x,y) to fractional grid indices, accounting for the
//      staggered half-cell offset of each field.
//   2. Split into integer cell index + fractional weight [0,1).
//   3. Clamp the integer index so we never read out of bounds.
//   4. Bilinear blend over the four surrounding nodes.

varType SemiLagrangian::interpolateU(const varType x, const varType y) const {
  // u is at (i.dx, (j+0.5).dy) → subtract the j offset before flooring.
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
  // v is at ((i+0.5)·dx, j·dy) → subtract the i offset before flooring.
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
