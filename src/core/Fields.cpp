#include "Fields.hpp"
#include <cmath>

void Fields2D::Div() {
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      const varType dudx = (u.Get(i + 1, j) - u.Get(i, j)) / dx;
      const varType dvdy = (v.Get(i, j + 1) - v.Get(i, j)) / dy;
      div.Set(i, j, dudx + dvdy);
    }
  }
}

void Fields2D::VelocityNormCenterGrid() {
  // Interpolate u and v from their staggered positions to cell centres, then
  // store the magnitude. The loop stops at nx-1 / ny-1 because the
  // cell-centre sample point (i + 0.5)*dx requires one ghost layer.
  for (int j = 0; j < ny - 1; j++) {
    for (int i = 0; i < nx - 1; i++) {
      const varType x = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
      const varType y = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;

      const varType uCenter = u.Interpolate(x, y, dx, dy, 0);
      const varType vCenter = v.Interpolate(x, y, dx, dy, 1);

      normVelocity.Set(i, j, std::sqrt(uCenter * uCenter + vCenter * vCenter));
    }
  }
}

void Fields2D::SolidCylinder(int cx, int cy, int r) {
  const int r2 = r * r;
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      const int ddx = i - cx;
      const int ddy = j - cy;
      if (ddx * ddx + ddy * ddy <= r2)
        SetLabel(i, j, SOLID);
    }
  }
}

void Fields2D::SolidBorders() {
  // Bottom and top rows — i is inner (fast index)
  for (int i = 0; i < nx; i++) {
    SetLabel(i, 0,      SOLID);
    SetLabel(i, ny - 1, SOLID);
  }
  // Left and right columns — j is inner (each call touches a different row,
  // stride is already 1 element in row-major so this is fine as-is)
  for (int j = 0; j < ny; j++) {
    SetLabel(0,      j, SOLID);
    SetLabel(nx - 1, j, SOLID);
  }
}
