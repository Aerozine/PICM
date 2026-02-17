#pragma once
#include "Grid2D.hpp"
#include <cstdint> // for uint8_t

class Fields2D {
public:
  enum CellType : uint8_t { FLUID = 0, SOLID = 1 };

  int nx, ny;
  varType density, dt, dx, dy;
  Grid2D u, v, p, div, rot, normVelocity;

  Fields2D(int nx, int ny, varType density, varType dt, varType dx,
           varType dy)
      : nx(nx), ny(ny), density(density), dt(dt), dx(dx), dy(dy),
        u(nx + 1, ny), v(nx, ny + 1), p(nx, ny), div(nx, ny),
        rot(nx, ny), normVelocity(nx, ny),  labels(nx * ny, FLUID) {}

  varType usolid = 0.0; // first try (fixed borders)
                        // next step -> moving boundaries ?

  CellType Label(int i, int j) const {
    return static_cast<CellType>(labels[idx(i, j)]);
  }

  void SetLabel(int i, int j, CellType t) {
    labels[idx(i, j)] = static_cast<uint8_t>(t);
  }

  void Div();
  void InitPotentialGradient();
  void InitTaylorGreen(const varType amplitude);
  void VelocityNormCenterGrid();
  void SolidCylinder(int cx, int cy, int r);
  void SolidBorders();

private:
  std::vector<uint8_t> labels; // not a Grid2D
  int idx(int i, int j) const { return ny * i + j; }
};
