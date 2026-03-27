#pragma once
#include "Grid2D.hpp"
#include <cstdint>
#include <string>
#include <vector>
#define FIELD_USOLID 0.0
class Fields2D {
public:
  /// @brief Possible states for a grid cell.

  enum CellType : uint16_t {
    FLUID = 0,      ///< Active fluid cell, participates in the pressure solve.
    SOLID = 1 << 0, ///< Solid (obstacle / wall) cell, velocity is fixed.

    BC_U =
        1
        << 1, ///< Boundary condition, horizontal velocity on the left is fixed.
    BC_V =
        1 << 2, ///< Boundary condition, vertical velocity underneath is fixed.
    BC_P = 1 << 3, ///< Boundary condition, pressure is fixed on cell center.
    BC_S = 1 << 4, ///< Boundary condition for smoke map.

    IC_U =
        1
        << 5, ///< Initial condition, horizontal velocity on the left is fixed.
    IC_V =
        1 << 6, ///< Initial condition, vertical velocity underneath is fixed.
    IC_P = 1 << 7, ///< Initial condition, pressure is fixed on cell center.
    IC_S = 1 << 8  ///< Initial condition for smoke map.
  };

  int nx,ny;
  varType density; ///< Fluid density (kg/m³).
  varType dt,dx,dy;
  Grid2D u,v,p,div;
  Grid2D normVelocity;
  Grid2D smokeMap;


    Grid2D u_sum, u_weight;
    Grid2D v_sum, v_weight;

    Grid2D countAliveParticles;
  // TODO improve string method to sth else and maybe particularize ?
  Fields2D(int nx, int ny, varType density, varType dt, varType dx, varType dy,
           std::string method)
      : nx(nx), ny(ny), density(density), dt(dt), dx(dx), dy(dy), u(nx + 1, ny),
        v(nx, ny + 1), p(nx, ny), div(nx, ny), normVelocity(nx - 1, ny - 1),
        smokeMap(nx - 1, ny - 1),
        u_sum(method == "PIC" ? nx + 1 : 0, method == "PIC" ? ny : 0),
        u_weight(method == "PIC" ? nx + 1 : 0, method == "PIC" ? ny : 0),
        v_sum(method == "PIC" ? nx : 0, method == "PIC" ? ny + 1 : 0),
        v_weight(method == "PIC" ? nx : 0, method == "PIC" ? ny + 1 : 0),
        countAliveParticles(method == "PIC" ? nx : 0, method == "PIC" ? ny : 0),
        labels(static_cast<std::size_t>(nx) * ny, FLUID) {}

  /// @brief Return the cell type (FLUID or SOLID) of cell (i, j).
  [[nodiscard]] CellType Label(int i, int j) const {
  return static_cast<CellType>(labels[idx(i, j)]);
  }

  /// @brief Set the cell type of cell (i, j).
  void SetLabel(int i, int j, CellType t) {
    labels[idx(i, j)] |= static_cast<uint16_t>(t);
  }

  ///@brief Compute the divergence to div
  void Div();

  /// @brief interpolate the velocity magnitude per cell to @c normVelocity
  void VelocityNormCenterGrid();

private:
  std::vector<uint16_t> labels; ///< Flat cell-type array, same layout as p.
  // TODO improve string method to sth else and maybe particularize ?

  [[nodiscard]] int idx(int i, int j) const { return nx * j + i; }
};
