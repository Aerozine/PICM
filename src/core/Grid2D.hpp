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
 * where values are written x-fastest, allowing the raw @c A buffer to be
 * passed directly to the writer without any transposition.
 *
 * All inner loops should therefore iterate over i in the innermost loop to
 * maximise cache locality.
 */
class Grid2D {
public:
  int nx; ///< Number of cells in the x-direction.
  int ny; ///< Number of cells in the y-direction.

  std::vector<varType> A; ///< Flat cell data, row-major: A[nx*j + i].

  /**
   * @brief Construct a zero-initialised grid of size @p nx × @p ny.
   */
  Grid2D(int nx, int ny) : nx(nx), ny(ny), A(nx * ny, varType{0}) {}

  /// @brief Read the value at cell (i, j).
  [[nodiscard]] varType Get(int i, int j) const { return A[nx * j + i]; }

  /// @brief Write a value into cell (i, j).
  void Set(int i, int j, varType val) { A[nx * j + i] = val; }

  /// @brief Return true if (i, j) is within grid bounds.
  [[nodiscard]] bool InBounds(int i, int j) const {
    return i >= 0 && i < nx && j >= 0 && j < ny;
  }

  /**
   * @brief Bilinearly interpolate this grid at physical position (x, y).
   *
   * Accounts for the staggered half-cell offset of each field type:
   * - @p field == 0 (u): nodes at (i·dx, (j+0.5)·dy) → j_real -= 0.5
   * - @p field == 1 (v): nodes at ((i+0.5)·dx, j·dy) → i_real -= 0.5
   * - Any other value: cell-centred, no offset.
   *
   * Indices are clamped so the 2×2 stencil always stays in bounds.
   *
   * @param x     Physical x-coordinate.
   * @param y     Physical y-coordinate.
   * @param dx    Cell width  in x (m).
   * @param dy    Cell height in y (m).
   * @param field Stagger type: 0 = u-face, 1 = v-face, other = cell-centre.
   */
  [[nodiscard]] varType Interpolate(varType x, varType y, varType dx,
                                    varType dy, int field) const;
};
