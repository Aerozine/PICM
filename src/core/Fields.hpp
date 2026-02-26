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
 * | Field         | Size           | Location                    |
 * |---------------|----------------|-----------------------------|
 * | @c u          | (nx+1) × ny    | x-face centres              |
 * | @c v          | nx × (ny+1)    | y-face centres              |
 * | @c p          | nx × ny        | cell centres                |
 * | @c div        | nx × ny        | cell centres (diagnostic)   |
 * | @c normVelocity | (nx-1) × (ny-1)      | cell centres (diagnostic)   |
 * | @c smokeMap | (nx-1) × (ny-1)      | cell centres (diagnostic)   |
 *
 * Cell labels (FLUID / SOLID) are stored in a separate flat array and
 * accessed via @c Label() / @c SetLabel().
 */
class Fields2D {
public:
  /// @brief Possible states for a grid cell.
  enum CellType : uint8_t {
    FLUID = 0, ///< Active fluid cell, participates in the pressure solve.
    SOLID = 1  ///< Solid (obstacle / wall) cell, velocity is fixed.
  };

  int nx;          ///< Number of pressure cells in x.
  int ny;          ///< Number of pressure cells in y.
  varType density; ///< Fluid density.
  varType dt;      ///< Time-step size.
  varType dx;      ///< Cell width  in x.
  varType dy;      ///< Cell height in y.

  Grid2D u;   ///< x-velocity, staggered: (nx+1) × ny.
  Grid2D v;   ///< y-velocity, staggered: nx × (ny+1).
  Grid2D p;   ///< Pressure,   cell-centred: nx × ny.
  Grid2D div; ///< Velocity divergence @f$\nabla \cdot \mathbf{u}$@f$
              ///< (diagnostic): @f$n_x \times n_y$@f$.
  Grid2D
      normVelocity; ///< |u| interpolated to cell centres (diagnostic): nx × ny.
  Grid2D smokeMap; ///< smoke matter in each cell centres

  /// Velocity imposed on SOLID cells (0 = no-slip). Reserved for moving
  /// boundaries in future work.
  varType usolid = REAL_LITERAL(0.0);

  /**
   * @brief Construct all fields and zero-initialise them.
   * @param nx      Number of pressure cells in x.
   * @param ny      Number of pressure cells in y.
   * @param density Fluid density.
   * @param dt      Time-step size.
   * @param dx      Cell width  in x.
   * @param dy      Cell height in y.
   */
  Fields2D(int nx, int ny, varType density, varType dt, varType dx, varType dy)
      : nx(nx), ny(ny), density(density), dt(dt), dx(dx), dy(dy), u(nx + 1, ny),
        v(nx, ny + 1), p(nx, ny), div(nx, ny), normVelocity(nx - 1, ny - 1),
        smokeMap(nx - 1, ny - 1),
        labels(static_cast<std::size_t>(nx) * ny, FLUID) {}

  // Cell label accessors
  /**
   * @brief Return the cell type (FLUID or SOLID) of cell (i, j).
   * @param i Column index [0, nx).
   * @param j Row    index [0, ny).
   */
  [[nodiscard]] CellType Label(int i, int j) const {
    return static_cast<CellType>(labels[idx(i, j)]);
  }

  /**
   * @brief Set the cell type of cell (i, j).
   * @param i Column index [0, nx).
   * @param j Row    index [0, ny).
   * @param t New cell type.
   */
  void SetLabel(int i, int j, CellType t) {
    labels[idx(i, j)] = static_cast<uint8_t>(t);
  }

  // Field update methods
  /**
   * @brief Compute the discrete divergence @f$\nabla \cdot \mathbf{u}$@f$ into
   * @c div.
   *
   * Uses first-order finite differences on the staggered grid:
   * @f[
   *   \text{div}(i,j) = \frac{u(i+1,j) - u(i,j)}{\Delta x}
   *                   + \frac{v(i,j+1) - v(i,j)}{\Delta y}
   * @f]
   */
  void Div();

  /**
   * @brief Interpolate the velocity magnitude |u| to cell centres and store
   *        the result in @c normVelocity.
   */
  void VelocityNormCenterGrid();

  // Geometry helpers

  /**
   * @brief Mark cells inside a disc as SOLID.
   * @param cx Centre x-index (cells).
   * @param cy Centre y-index (cells).
   * @param r  Radius in cells.
   */
  void SolidCylinder(int cx, int cy, int r);

  /**
   * @brief Mark the four border rows/columns as SOLID (no-slip walls).
   */
  void SolidBorders();

private:
  std::vector<uint8_t> labels; ///< Flat cell-type array, same layout as p.

  /// @brief Flat index into @c labels (column-major, matching Grid2D).
  [[nodiscard]] int idx(int i, int j) const { return ny * i + j; }
};
