#pragma once
#include "Grid2D.hpp"
#include <cstdint> // for uint8_t

class Fields2D {
public:
  enum CellType : uint8_t { FLUID = 0, SOLID = 1 };

  int nx, ny;
  varType density, dt, dx, dy;
  Grid2D u, v, p, div, rot;

  Fields2D(int nx, int ny, varType density, varType dt, varType dx,
           varType dy)
      : nx(nx), ny(ny), density(density), dt(dt), dx(dx), dy(dy),
        u(nx - 1, ny), // !
        v(nx, ny - 1), p(nx - 1, ny - 1), div(nx - 1, ny - 1),
        rot(nx - 1, ny - 1), labels((nx - 1) * (ny - 1), FLUID) {}

  varType usolid = 0.0; // first try (fixed borders)
                        // next step -> moving boundaries ?
// nodiscard for clang-tidy not usefull to understand
  [[nodiscard]] CellType Label(int i, int j) const {
    return static_cast<CellType>(labels[idx(i, j)]);
  }

  void SetLabel(int i, int j, CellType t) {
    labels[idx(i, j)] = static_cast<uint8_t>(t);
  }

  void Div();
  void InitPotentialGradient(varType amplitude, int kx, int ky);
  void InitTaylorGreen(varType amplitude);
  // Grid2D VelocityNormCenterGrid();

private:
  std::vector<uint8_t> labels; // not a Grid2D
  // FIX: Use (nx-1) as stride since labels is (nx-1) x (ny-1)
  [[nodiscard]] int idx(int i, int j) const { return (nx - 1) * j + i; }
};
