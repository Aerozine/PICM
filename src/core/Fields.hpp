#pragma once
/**
 * @file Fields.hpp
 * @brief Physical fields for a 2-D incompressible simulation on a MAC grid.
 */
#include "Grid2D.hpp"
#include <cstdint>
#include <vector>

class Fields2D {
public:
  /// @brief Possible states for a grid cell.
  enum CellType : uint8_t {
    FLUID = 0, ///< Active fluid cell, participates in the pressure solve.
    SOLID = 1  ///< Solid (obstacle / wall) cell, velocity is fixed to usolid.
  };

  int nx,ny;
  varType density; ///< Fluid density (kg/m³).
  varType dt,dx,dy;
  Grid2D u,v,p,div;
  Grid2D normVelocity;
  Grid2D smokeMap;


  /**
   * @brief Construct all fields and zero-initialise them.
   */
  Fields2D(int nx, int ny, varType density, varType dt, varType dx, varType dy)
      : nx(nx), ny(ny), density(density), dt(dt), dx(dx), dy(dy), u(nx + 1, ny),
        v(nx, ny + 1), p(nx, ny), div(nx, ny), normVelocity(nx, ny),
        smokeMap(nx, ny), labels(static_cast<std::size_t>(nx) * ny, FLUID) {}

  /// @brief Return the cell type (FLUID or SOLID) of cell (i, j).
  [[nodiscard]] CellType Label(int i, int j) const {
  return static_cast<CellType>(labels[idx(i, j)]);
  }

  /// @brief Set the cell type of cell (i, j).
  void SetLabel(int i, int j, CellType t) {
    labels[idx(i, j)] = static_cast<uint8_t>(t);
  }

  ///@brief Compute the divergence to div
  void Div();

  /// @brief interpolate the velocity magnitude per cell to @c normVelocity
  void VelocityNormCenterGrid();

private:
  std::vector<uint8_t> labels; ///< Flat cell-type array, same layout as p.

  [[nodiscard]] int idx(int i, int j) const { return nx * j + i; }
};
