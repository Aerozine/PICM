#include "Fields.hpp"

#include <cassert>
#include <cmath>

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

      const varType uc = u.interpolate(x, y, dx, dy, 0);
      const varType vc = v.interpolate(x, y, dx, dy, 1);

      normVelocity.Set(i, j, std::sqrt(uc * uc + vc * vc));
    }
}
