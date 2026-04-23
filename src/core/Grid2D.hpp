#pragma once
#include <cstdlib>
#include <cstring>

#include "Precision.hpp"
#include <cassert>
#include <cmath>
#include <vector>
///@brief  **ROW-MAJOR 2d grid** @c A[nx*j+i]
class Grid2D {
  //@todo do a sparsified version and an eigen version if eigen defined

public:
  int nx;
  int ny;
  Grid2D()=default;

  // calloc is just malloc but init at 0
  Grid2D(int nx_, int ny_)
      : nx(nx_), ny(ny_),
        A(static_cast<varType *>(
            std::calloc(static_cast<std::size_t>(nx_) * ny_, sizeof(varType))))
  {}
  ~Grid2D() { std::free(A); }
  Grid2D &operator=(const Grid2D &o) {
    if (this == &o) return *this;
    if (nx * ny != o.nx * o.ny) {
      std::free(A);
      A = static_cast<varType *>(
          std::malloc(static_cast<std::size_t>(o.nx) * o.ny * sizeof(varType)));
    }
    nx = o.nx; ny = o.ny;
    std::memcpy(A, o.A, static_cast<std::size_t>(nx) * ny * sizeof(varType));
    return *this;
  }
  Grid2D(Grid2D &&o) noexcept : nx(o.nx), ny(o.ny), A(o.A) {
    o.nx = 0; o.ny = 0; o.A = nullptr;
  }
  // speed up due to AVX/SSE
  void reset() noexcept {
    std::memset(A, 0, static_cast<std::size_t>(nx) * ny * sizeof(varType));
  }
  Grid2D(const Grid2D& o)
    : nx(o.nx), ny(o.ny),
      A(static_cast<varType*>(std::malloc(
            static_cast<std::size_t>(o.nx) * o.ny * sizeof(varType))))
  {
    std::memcpy(A, o.A, static_cast<std::size_t>(nx) * ny * sizeof(varType));
  }
  Grid2D &operator=(Grid2D &&o) noexcept {
    if (this == &o) return *this;
    std::free(A);
    nx = o.nx; ny = o.ny; A = o.A;
    o.nx = 0; o.ny = 0; o.A = nullptr;
    return *this;
  }
  /// @brief Read the value at cell (i, j).
  [[nodiscard]] varType Get(const int i, const int j) const {
    return A[nx * j + i];
  }

  /// @brief Write a value into cell (i, j).
  void Set(const int i, const int j, const varType val) { A[nx * j + i] = val; }

  /**
   * @brief Bilinearly interpolation
   * @param x     Physical x-coordinate.
   * @param y     Physical y-coordinate.
   * @param dx    Cell width  in x.
   * @param dy    Cell height in y.
   * @param field Stagger type: 0 = u-face, 1 = v-face, other = cell-centre.
   */
  [[nodiscard]] varType interpolate(varType x, varType y, varType dx,
                                    varType dy, int field) const;
//private:

  varType *A = nullptr;
};
