#pragma once
/**
 * @file Fields.hpp
 * @brief Physical fields for a 2-D incompressible simulation on a MAC grid.
 */
#include "Grid2D.hpp"
#include <cstdint>
#include <string>
#include <vector>

class Fields2D {
public:
  /// @brief Possible states for a grid cell (bitmask).
  enum CellType : uint16_t {
    FLUID = 0,      ///< Active fluid cell, participates in the pressure solve.
    SOLID = 1 << 0, ///< Solid (obstacle / wall) cell, velocity is fixed.

    BC_U = 1 << 1, ///< Boundary condition: horizontal velocity on the left is fixed.
    BC_V = 1 << 2, ///< Boundary condition: vertical velocity underneath is fixed.
    BC_P = 1 << 3, ///< Boundary condition: pressure is fixed on cell center.
    BC_S = 1 << 4, ///< Boundary condition for smoke map.

    IC_U = 1 << 5, ///< Initial condition: horizontal velocity on the left is fixed.
    IC_V = 1 << 6, ///< Initial condition: vertical velocity underneath is fixed.
    IC_P = 1 << 7, ///< Initial condition: pressure is fixed on cell center.
    IC_S = 1 << 8  ///< Initial condition for smoke map.
  };

  int nx, ny;
  varType density; ///< Fluid density (kg/m3).
  varType dt, dx, dy;

  Grid2D u;            ///< x-velocity, staggered: (nx+1) x ny.
  Grid2D v;            ///< y-velocity, staggered: nx x (ny+1).
  Grid2D p;            ///< Pressure, cell-centred: nx x ny.
  Grid2D div;          ///< Velocity divergence (diagnostic): nx x ny.
  Grid2D normVelocity; ///< |u| interpolated to cell centres (diagnostic): nx x ny.
  Grid2D smokeMap;     ///< Smoke density at cell centres: nx x ny.

  /// PIC-only P2G accumulation buffers (zero-sized for SemiLagrangian).
  Grid2D u_sum, u_weight;
  Grid2D v_sum, v_weight;
  Grid2D countAliveParticles;

  /// Velocity imposed on SOLID cells (0 = no-slip).
  varType usolid = REAL_LITERAL(0.0);

  /**
   * @brief Construct all fields and zero-initialise them.
   * @param method "PIC" allocates the P2G accumulation buffers; anything else
   *               leaves them zero-sized.
   */
  Fields2D(int nx, int ny, varType density, varType dt, varType dx, varType dy,
           std::string method)
      : nx(nx), ny(ny), density(density), dt(dt), dx(dx), dy(dy),
        u(nx + 1, ny), v(nx, ny + 1), p(nx, ny), div(nx, ny),
        normVelocity(nx, ny), smokeMap(nx, ny),
        u_sum    (method == "PIC" ? nx + 1 : 0, method == "PIC" ? ny     : 0),
        u_weight (method == "PIC" ? nx + 1 : 0, method == "PIC" ? ny     : 0),
        v_sum    (method == "PIC" ? nx     : 0, method == "PIC" ? ny + 1 : 0),
        v_weight (method == "PIC" ? nx     : 0, method == "PIC" ? ny + 1 : 0),
        countAliveParticles(method == "PIC" ? nx : 0, method == "PIC" ? ny : 0),
        labels(static_cast<std::size_t>(nx) * ny, FLUID) {}

  /// @brief Return the cell label bitmask of cell (i, j).
  [[nodiscard]] uint16_t Label(int i, int j) const {
    return labels[idx(i, j)];
  }

  /// @brief OR a cell type flag into cell (i, j).
  void SetLabel(int i, int j, CellType t) {
    labels[idx(i, j)] |= static_cast<uint16_t>(t);
  }

  /// @brief Compute the velocity divergence into @c div.
  void Div();

  /// @brief Interpolate the velocity magnitude per cell into @c normVelocity.
  void VelocityNormCenterGrid();

private:
  std::vector<uint16_t> labels; ///< Flat cell-type bitmask, same layout as p.

  [[nodiscard]] int idx(int i, int j) const { return nx * j + i; }
};
