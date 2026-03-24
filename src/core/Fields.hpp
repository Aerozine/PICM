#pragma once
#include "Grid2D.hpp"
#include <cstdint>
#include <vector>

/**
 * @file Fields.hpp
 * @brief Physical fields for a 2-D incompressible simulation on a MAC grid.
 */

/**
 * @brief All physical fields for a 2-D incompressible Navier-Stokes solver
 *        on a staggered (MAC / Marker-And-Cell) grid.
 *
 * ### Grid layout
 * | Field           | Size        | Location                         |
 * |-----------------|-------------|----------------------------------|
 * | @c u            | (nx+1) × ny | x-face centres                   |
 * | @c v            | nx × (ny+1) | y-face centres                   |
 * | @c p            | nx × ny     | cell centres                     |
 * | @c div          | nx × ny     | cell centres (diagnostic)        |
 * | @c normVelocity | nx × ny     | cell centres (diagnostic)        |
 * | @c smokeMap     | nx × ny     | cell centres                     |
 *
 * Cell labels (FLUID / SOLID) are stored in a separate flat array and
 * accessed via @c Label() / @c SetLabel().
 */
class Fields2D {
public:
  /// @brief Possible states for a grid cell.
  enum CellType : uint8_t {
    FLUID = 0, ///< Active fluid cell, participates in the pressure solve.
    SOLID = 1  ///< Solid (obstacle / wall) cell, velocity is fixed to usolid.
  };

  int nx;          ///< Number of pressure cells in x.
  int ny;          ///< Number of pressure cells in y.
  varType density; ///< Fluid density (kg/m³).
  varType dt;      ///< Time-step size (s).
  varType dx;      ///< Cell width  in x (m).
  varType dy;      ///< Cell height in y (m).

  Grid2D u;            ///< x-velocity, staggered: (nx+1) × ny.
  Grid2D v;            ///< y-velocity, staggered: nx × (ny+1).
  Grid2D p;            ///< Pressure, cell-centred: nx × ny.
  Grid2D div;          ///< Velocity divergence ∇·u (diagnostic): nx × ny.
  Grid2D normVelocity; ///< |u| at cell centres (diagnostic): nx × ny.
  Grid2D smokeMap;     ///< Smoke concentration, cell-centred: nx × ny.

  /// Velocity imposed on SOLID faces (0 = no-slip).
  varType usolid = REAL_LITERAL(0.0);

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

  /**
   * @brief Compute the discrete divergence ∇·u into @c div.
   *
   * \f$ \mathrm{div}(i,j)
   *   = \frac{u(i+1,j) - u(i,j)}{\Delta x}
   *   + \frac{v(i,j+1) - v(i,j)}{\Delta y} \f$
   */
  void Div();

  /**
   * @brief Interpolate the velocity magnitude |u| to cell centres and store
   *        the result in @c normVelocity.
   */
  void VelocityNormCenterGrid();

private:
  std::vector<uint8_t> labels; ///< Flat cell-type array, same layout as p.

  [[nodiscard]] int idx(int i, int j) const { return nx * j + i; }
};
