#pragma once
#include "Grid2D.hpp"
#include <cstdint> // for uint8_t

class Fields2D {
public:
  enum CellType : uint8_t { FLUID = 0, SOLID = 1 };

  size_t nx, ny;
  Real density, dt, dx, dy;
  Grid2D u, v, p, div, rot;

  Fields2D(size_t nx, size_t ny, Real density, Real dt, Real dx, Real dy)
      : nx(nx), ny(ny), density(density), dt(dt), dx(dx), dy(dy), u(nx, ny - 1),
        v(nx - 1, ny), p(nx - 1, ny - 1), div(nx - 1, ny - 1),
        rot(nx - 1, ny - 1), labels((nx - 1) * (ny - 1), FLUID) {}

  Real usolid = 0.0; // first try (fixed borders)
                     // next step -> moving boundaries ?

  CellType Label(size_t i, size_t j) const {
    return static_cast<CellType>(labels[idx(i, j)]);
  }

  void SetLabel(size_t i, size_t j, CellType t) {
    labels[idx(i, j)] = static_cast<uint8_t>(t);
  }

  void Div();
  void InitRandomVelocities();
  void InitPotentialGradient(Real amplitude, int kx, int ky);
  Grid2D VelocityNormCenterGrid();

private:
  std::vector<uint8_t> labels; // not a Grid2D
  // FIXED: Use (nx-1) as stride since labels has dimensions (nx-1) x (ny-1)
  size_t idx(size_t i, size_t j) const { return (nx - 1) * j + i; }
};
