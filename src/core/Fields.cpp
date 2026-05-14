#include "Fields.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>


varType bilerp(const varType f00, const varType f10,
                             const varType f01, const varType f11,
                             const varType fx, const varType fy) {
  return (REAL_LITERAL(1.0) - fy) *
             ((REAL_LITERAL(1.0) - fx) * f00 + fx * f10) +
         fy * ((REAL_LITERAL(1.0) - fx) * f01 + fx * f11);
}

varType sampleUFreeSlip(const Fields2D &fields,
                                      const Grid2D &grid, int i, int j) {
  i = std::clamp(i, 0, grid.nx - 1);
  j = std::clamp(j, 0, grid.ny - 1);

  const bool leftSolid = IS_SOLID(fields.Label(i, j + 1));
  const bool rightSolid = IS_SOLID(fields.Label(i + 1, j + 1));
  if (leftSolid != rightSolid)
    return REAL_LITERAL(0.0);
  if (!leftSolid)
    return grid.Get(i, j);

  if (j > 0 && !IS_SOLID(fields.Label(i, j)) &&
      !IS_SOLID(fields.Label(i + 1, j)))
    return grid.Get(i, j - 1);
  if (j + 1 < grid.ny && !IS_SOLID(fields.Label(i, j + 2)) &&
      !IS_SOLID(fields.Label(i + 1, j + 2)))
    return grid.Get(i, j + 1);
  return REAL_LITERAL(0.0);
}

varType sampleVFreeSlip(const Fields2D &fields,
                                      const Grid2D &grid, int i, int j) {
  i = std::clamp(i, 0, grid.nx - 1);
  j = std::clamp(j, 0, grid.ny - 1);

  const bool bottomSolid = IS_SOLID(fields.Label(i + 1, j));
  const bool topSolid = IS_SOLID(fields.Label(i + 1, j + 1));
  if (bottomSolid != topSolid)
    return REAL_LITERAL(0.0);
  if (!bottomSolid)
    return grid.Get(i, j);

  if (i > 0 && !IS_SOLID(fields.Label(i, j)) &&
      !IS_SOLID(fields.Label(i, j + 1)))
    return grid.Get(i - 1, j);
  if (i + 1 < grid.nx && !IS_SOLID(fields.Label(i + 2, j)) &&
      !IS_SOLID(fields.Label(i + 2, j + 1)))
    return grid.Get(i + 1, j);
  return REAL_LITERAL(0.0);
}


varType Fields2D::interpolateU(varType x, varType y) const {
  return interpolateU(u, x, y);
}

varType Fields2D::interpolateV(varType x, varType y) const {
  return interpolateV(v, x, y);
}

varType Fields2D::interpolateU(const Grid2D &grid, const varType x,
                               const varType y) const {
  assert(grid.nx == nx + 1);
  assert(grid.ny == ny);

  const varType iReal =
      std::clamp(x / dx, REAL_LITERAL(0.0), static_cast<varType>(grid.nx - 1));
  const varType jReal =
      std::clamp(y / dy - REAL_LITERAL(0.5), REAL_LITERAL(0.0),
                 static_cast<varType>(grid.ny - 1));

  const int i0 = static_cast<int>(std::floor(iReal));
  const int j0 = static_cast<int>(std::floor(jReal));
  const int i1 = std::min(i0 + 1, grid.nx - 1);
  const int j1 = std::min(j0 + 1, grid.ny - 1);

  const varType value =
      bilerp(sampleUFreeSlip(*this, grid, i0, j0),
             sampleUFreeSlip(*this, grid, i1, j0),
             sampleUFreeSlip(*this, grid, i0, j1),
             sampleUFreeSlip(*this, grid, i1, j1), iReal - i0, jReal - j0);
  assert(std::isfinite(value));
  return value;
}

varType Fields2D::interpolateV(const Grid2D &grid, const varType x,
                               const varType y) const {
  assert(grid.nx == nx);
  assert(grid.ny == ny + 1);

  const varType iReal =
      std::clamp(x / dx - REAL_LITERAL(0.5), REAL_LITERAL(0.0),
                 static_cast<varType>(grid.nx - 1));
  const varType jReal =
      std::clamp(y / dy, REAL_LITERAL(0.0), static_cast<varType>(grid.ny - 1));

  const int i0 = static_cast<int>(std::floor(iReal));
  const int j0 = static_cast<int>(std::floor(jReal));
  const int i1 = std::min(i0 + 1, grid.nx - 1);
  const int j1 = std::min(j0 + 1, grid.ny - 1);

  const varType value =
      bilerp(sampleVFreeSlip(*this, grid, i0, j0),
             sampleVFreeSlip(*this, grid, i1, j0),
             sampleVFreeSlip(*this, grid, i0, j1),
             sampleVFreeSlip(*this, grid, i1, j1), iReal - i0, jReal - j0);
  assert(std::isfinite(value));
  return value;
}

void Fields2D::UpdateDivNorm() {
    Div();
    VelocityNormCenterGrid();
}

void Fields2D::Div() {
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      if (IS_AIR(Label(i + 1, j + 1))) {
        div.Set(i, j, 0.0);
        continue;
      }

      const varType dudx = (u.Get(i + 1, j) - u.Get(i, j)) / dx;
      const varType dvdy = (v.Get(i, j + 1) - v.Get(i, j)) / dy;
      assert(std::isfinite(dudx));
      assert(std::isfinite(dvdy));
      div.Set(i, j, dudx + dvdy);
    }
}

void Fields2D::VelocityNormCenterGrid() {
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const varType x = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
      const varType y = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;

      const varType uc = interpolateU(x, y);
      const varType vc = interpolateV(x, y);

      normVelocity.Set(i, j, std::sqrt(uc * uc + vc * vc));
    }
}


void Fields2D::VorticityCenterGrid() {
  const varType xMax = static_cast<varType>(nx) * dx;
  const varType yMax = static_cast<varType>(ny) * dy;

  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const varType x = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
      const varType y = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;

      const varType xLeft =
          std::max(REAL_LITERAL(0.0), x - dx * REAL_LITERAL(0.5));
      const varType xRight = std::min(xMax, x + dx * REAL_LITERAL(0.5));
      const varType yBottom =
          std::max(REAL_LITERAL(0.0), y - dy * REAL_LITERAL(0.5));
      const varType yTop = std::min(yMax, y + dy * REAL_LITERAL(0.5));

      const varType dvdx =
          (interpolateV(xRight, y) - interpolateV(xLeft, y)) /
          std::max(xRight - xLeft, REAL_EPSILON);
      const varType dudy =
          (interpolateU(x, yTop) - interpolateU(x, yBottom)) /
          std::max(yTop - yBottom, REAL_EPSILON);

      vorticity.Set(i, j, dvdx - dudy);
    }
}