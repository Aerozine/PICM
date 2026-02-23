#include "Fields.hpp"
#include <cmath>

void Fields2D::Div() {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
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
  for (int i = 0; i < nx - 1; i++) {
    for (int j = 0; j < ny - 1; j++) {
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
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      const int ddx = i - cx;
      const int ddy = j - cy;
      if (ddx * ddx + ddy * ddy <= r2)
        SetLabel(i, j, SOLID);
    }
  }
}

void Fields2D::SolidBorders() {
  // Bottom and top rows
  for (int i = 0; i < nx; i++) {
    SetLabel(i, 0, SOLID);
    SetLabel(i, ny - 1, SOLID);
  }
  // Left and right columns
  for (int j = 0; j < ny; j++) {
    SetLabel(0, j, SOLID);
    SetLabel(nx - 1, j, SOLID);
  }
}
