#include "Fields.hpp"
#include <cmath>

void Fields2D::Div() {
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const varType dudx = (u.Get(i + 1, j) - u.Get(i, j)) / dx;
      const varType dvdy = (v.Get(i, j + 1) - v.Get(i, j)) / dy;
      div.Set(i, j, dudx + dvdy);
    }
}

void Fields2D::VelocityNormCenterGrid() {
  // normVelocity is nx × ny (cell-centred).
  // For cell (i, j), the centre is at ((i+0.5)*dx, (j+0.5)*dy).
  // We interpolate u and v from their staggered positions to that point.
  // The last cell column/row is included - Grid2D::Interpolate() clamps the
  // stencil, so there is no out-of-bounds access.

  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const varType x = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
      const varType y = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;

      const varType uc = u.Interpolate(x, y, dx, dy, 0);
      const varType vc = v.Interpolate(x, y, dx, dy, 1);

      normVelocity.Set(i, j, std::sqrt(uc * uc + vc * vc));
    }
}
