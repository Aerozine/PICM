#pragma once
#include "Precision.hpp"
#include <vector>

/**
 * @file Grid2D.hpp
 * @brief 2D scalar grid on a structured Cartesian mesh.
 */

/**
 * @brief A flat, heap-allocated 2D scalar grid.
 *
 * Data is stored in **row-major** order: element (i, j) lives at
 * @c A[nx * j + i], so the i-index (x-direction) is the fast index.
 *
 * This layout matches the VTK ImageData convention for appended binary data,
 * where values are written x-fastest (for z … for y … for x), allowing the
 * raw @c A buffer to be passed directly to the writer without any
 * transposition.
 *
 * All inner loops should therefore iterate over i in the innermost loop to
 * maximise cache locality.
 *
 * Grid dimensions are runtime values (read from a JSON config), so the
 * storage uses @c std::vector which is equivalent to a raw heap allocation
 * but provides automatic memory management and bounds-checking in debug builds.
 */
class Grid2D {
public:
  int nx; ///< Number of cells in the x-direction.
  int ny; ///< Number of cells in the y-direction.

  std::vector<varType> A; ///< Flat cell data, row-major: A[nx*j + i].

  /**
   * @brief Construct a zero-initialised grid of size @p nx × @p ny.
   * @param nx Number of cells in x.
   * @param ny Number of cells in y.
   */
  Grid2D(int nx, int ny) : nx(nx), ny(ny), A(nx * ny, varType{0}) {}

  /**
   * @brief Read the scalar value stored at cell (i, j).
   * @param i Column index (x), must be in [0, nx).
   * @param j Row    index (y), must be in [0, ny).
   * @return  Value at (i, j).
   */
  [[nodiscard]] varType Get(int i, int j) const { return A[nx * j + i]; }

  /**
   * @brief Write a scalar value into cell (i, j).
   * @param i   Column index (x), must be in [0, nx).
   * @param j   Row    index (y), must be in [0, ny).
   * @param val Value to store.
   */
  void Set(int i, int j, varType val) { A[nx * j + i] = val; }

  /**
   * @brief Check whether indices (i, j) lie inside the grid.
   * @param i Column index.
   * @param j Row    index.
   * @return  @c true if both indices are in bounds.
   */
  [[nodiscard]] bool InBounds(int i, int j) const {
    return i >= 0 && i < nx && j >= 0 && j < ny;
  }

  /**
   * @brief Bilinearly interpolate this grid at a physical position (x, y).
   *
   * Accounts for the staggered half-cell offset of each field type:
   * - @p field == 0 (u): node positions are (i·dx, (j+0.5)·dy) → subtract
   *   0.5 from the j fractional index.
   * - @p field == 1 (v): node positions are ((i+0.5)·dx, j·dy) → subtract
   *   0.5 from the i fractional index.
   * - Any other value: cell-centred, no offset applied.
   *
   * Indices are clamped so the four-node stencil always stays in bounds.
   *
   * @param x     Physical x-coordinate.
   * @param y     Physical y-coordinate.
   * @param dx    Cell width  in x (m).
   * @param dy    Cell height in y (m).
   * @param field Stagger type: 0 = u-face, 1 = v-face, other = cell-centre.
   * @return      Interpolated value.
   */
  [[nodiscard]] varType Interpolate(varType x, varType y, varType dx,
                                    varType dy, int field) const;
};
