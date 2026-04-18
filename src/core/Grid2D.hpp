#pragma once
#include "Precision.hpp"
#include <vector>
///@brief  **ROW-MAJOR 2d grid** @c A[nx*j+i]
class Grid2D {
public:
  int nx;
  int ny;
  /* for debuging purpose
  //  varType * A;
  //Grid2D(int nx, int ny) : nx(nx), ny(ny) {
  //  A = new varType[nx * ny];
  //}
  */
  std::vector<varType> A; ///< Flat cell data, row-major: A[nx*j + i].
  Grid2D(int nx, int ny) : nx(nx), ny(ny), A(nx * ny, varType{0}) {}

  /// @brief Read the value at cell (i, j).
  [[nodiscard]] varType Get(const int i, const int j) const {
    return A[nx * j + i];
  }

  /// @brief Write a value into cell (i, j).
  void Set(const int i, const int j, const varType val) { A[nx * j + i] = val; }

  /// @brief Return true if (i, j) is within grid bounds.
  [[nodiscard]] bool InBounds(int i, int j) const {
    return i >= 0 && i < nx && j >= 0 && j < ny;
  }

  /**
   * @brief Bilinearly interpolation
   * @param x     Physical x-coordinate.
   * @param y     Physical y-coordinate.
   * @param dx    Cell width  in x.
   * @param dy    Cell height in y.
   * @param field Stagger type: 0 = u-face, 1 = v-face, other = cell-centre.
   */
  [[nodiscard]] varType Interpolate(varType x, varType y, varType dx,
                                    varType dy, int field) const;
};
