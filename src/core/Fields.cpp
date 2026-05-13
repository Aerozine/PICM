#include "Fields.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace {


[[nodiscard]] varType sampleUFreeSlip(const Fields2D &fields,
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

[[nodiscard]] varType sampleVFreeSlip(const Fields2D &fields,
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

} // namespace

template <unsigned char field>
varType Fields2D::interpolate(const Grid2D &grid, varType x, varType y) const {
  varType i_real = x / dx;
  varType j_real = y / dy;

  if constexpr (field == 0) {
    // u-face: no x offset (face at i·dx), half-cell offset in y
    j_real -= REAL_LITERAL(0.5);
  } else if constexpr (field == 1) {
    // v-face: half-cell offset in x, no y offset (face at j·dy)
    i_real -= REAL_LITERAL(0.5);
  } else {
    // cell-centre: half-cell offset in both axes
    i_real -= REAL_LITERAL(0.5);
    j_real -= REAL_LITERAL(0.5);
  }

  const varType ic =
      std::clamp(i_real, REAL_LITERAL(0.0), static_cast<varType>(grid.nx - 1));
  const varType jc =
      std::clamp(j_real, REAL_LITERAL(0.0), static_cast<varType>(grid.ny - 1));

  const int i0 = static_cast<int>(std::floor(ic));
  const int j0 = static_cast<int>(std::floor(jc));
  const int i1 = std::min(i0 + 1, grid.nx - 1);
  const int j1 = std::min(j0 + 1, grid.ny - 1);

  const varType fx = ic - static_cast<varType>(i0);
  const varType fy = jc - static_cast<varType>(j0);

  varType f00, f10, f01, f11;
  if constexpr (field == 0) {
    f00 = sampleUFreeSlip(*this, grid, i0, j0);
    f10 = sampleUFreeSlip(*this, grid, i1, j0);
    f01 = sampleUFreeSlip(*this, grid, i0, j1);
    f11 = sampleUFreeSlip(*this, grid, i1, j1);
  } else if constexpr (field == 1) {
    f00 = sampleVFreeSlip(*this, grid, i0, j0);
    f10 = sampleVFreeSlip(*this, grid, i1, j0);
    f01 = sampleVFreeSlip(*this, grid, i0, j1);
    f11 = sampleVFreeSlip(*this, grid, i1, j1);
  } else {
    f00 = grid.Get(i0, j0);
    f10 = grid.Get(i1, j0);
    f01 = grid.Get(i0, j1);
    f11 = grid.Get(i1, j1);
  }

  const varType value =
      (REAL_LITERAL(1.0) - fy) *
          ((REAL_LITERAL(1.0) - fx) * f00 + fx * f10) +
      fy * ((REAL_LITERAL(1.0) - fx) * f01 + fx * f11);
  assert(std::isfinite(value));
  return value;
}

template varType Fields2D::interpolate<0>(const Grid2D &, varType, varType) const;
template varType Fields2D::interpolate<1>(const Grid2D &, varType, varType) const;
template varType Fields2D::interpolate<2>(const Grid2D &, varType, varType) const;

void Fields2D::UpdateDivNorm() {
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const varType x = (REAL_LITERAL(0.5) + static_cast<float>(i)) * dx;
      const varType y = (REAL_LITERAL(0.5) + static_cast<float>(j)) * dy;
      const varType uc = interpolateU(x, y);
      const varType vc = interpolateV(x, y);
      normVelocity.Set(i, j, std::sqrt(uc * uc + vc * vc));
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
