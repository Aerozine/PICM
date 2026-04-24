#include "Fields.hpp"

#include <cassert>
#include <cmath>
void Fields2D::UpdateDivNorm() {
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const varType x = (REAL_LITERAL(0.5) + static_cast<float>(i)) * dx;
      const varType y = (REAL_LITERAL(0.5) + static_cast<float>(j)) * dy;
      const varType uc = u.interpolate<0>(x, y, dx, dy);
      const varType vc = v.interpolate<1>(x, y, dx, dy);
      normVelocity.Set(i, j, std::sqrt(uc * uc + vc * vc));
      if (IS_AIR(Label(i + 1, j + 1))){
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
      if (IS_AIR(Label(i + 1, j + 1))){
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

      const varType uc = u.interpolate<0>(x, y, dx, dy);
      const varType vc = v.interpolate<1>(x, y, dx, dy);

      normVelocity.Set(i, j, std::sqrt(uc * uc + vc * vc));
    }
}
