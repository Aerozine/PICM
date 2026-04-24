#pragma once
#include <algorithm>
#include <cstdlib>
#include <cstring>

#include "Precision.hpp"
#include <cassert>
#include <cmath>
#include <vector>
///@brief  **ROW-MAJOR 2d grid** @c A[nx*j+i]
class Grid2D {
public:
  int nx;
  int ny;
  Grid2D() = default;

  // calloc is just malloc but init at 0
  Grid2D(int nx_, int ny_)
      : nx(nx_), ny(ny_),
        A(static_cast<varType *>(std::calloc(
            static_cast<std::size_t>(nx_) * ny_, sizeof(varType)))) {}
  ~Grid2D() { std::free(A); }
  Grid2D &operator=(const Grid2D &o) {
    if (this == &o)
      return *this;
    if (nx * ny != o.nx * o.ny) {
      std::free(A);
      A = static_cast<varType *>(
          std::malloc(static_cast<std::size_t>(o.nx) * o.ny * sizeof(varType)));
    }
    nx = o.nx;
    ny = o.ny;
    std::memcpy(A, o.A, static_cast<std::size_t>(nx) * ny * sizeof(varType));
    return *this;
  }
  Grid2D(Grid2D &&o) noexcept : nx(o.nx), ny(o.ny), A(o.A) {
    o.nx = 0;
    o.ny = 0;
    o.A = nullptr;
  }
  // speed up due to AVX/SSE
  void reset() noexcept {
    std::memset(A, 0, static_cast<std::size_t>(nx) * ny * sizeof(varType));
  }
  Grid2D(const Grid2D &o)
      : nx(o.nx), ny(o.ny),
        A(static_cast<varType *>(std::malloc(static_cast<std::size_t>(o.nx) *
                                             o.ny * sizeof(varType)))) {
    std::memcpy(A, o.A, static_cast<std::size_t>(nx) * ny * sizeof(varType));
  }
  [[nodiscard]] inline varType &operator()(int i, int j) noexcept {
    return A[nx * j + i];
  }
  [[nodiscard]] inline varType operator()(int i, int j) const noexcept {
    return A[nx * j + i];
  }
  Grid2D &operator=(Grid2D &&o) noexcept {
    if (this == &o)
      return *this;
    std::free(A);
    nx = o.nx;
    ny = o.ny;
    A = o.A;
    o.nx = 0;
    o.ny = 0;
    o.A = nullptr;
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
  // template is at compile time , so if we
  // put 0 in field , it simplify directly in compilation
  // same for 1 and 2
  template <__uint8_t field>
  varType interpolate(const varType x, const varType y, const varType dx,
                      const varType dy) const {
    varType i_real = x / dx;
    varType j_real = y / dy;

    if constexpr (field == 0)
      j_real -= REAL_LITERAL(0.5); // u-face: staggered in y
    if constexpr (field == 1)
      i_real -= REAL_LITERAL(0.5); // v-face: staggered in x
    // else no staggering

    const varType i_clamped =
        std::clamp(i_real, REAL_LITERAL(0.0), static_cast<varType>(nx - 1));
    const varType j_clamped =
        std::clamp(j_real, REAL_LITERAL(0.0), static_cast<varType>(ny - 1));

    const int i0 = static_cast<int>(std::floor(i_clamped));
    const int j0 = static_cast<int>(std::floor(j_clamped));

    const varType fx = i_clamped - static_cast<varType>(i0);
    const varType fy = j_clamped - static_cast<varType>(j0);

    const int i1 = std::min(i0 + 1, nx - 1);
    const int j1 = std::min(j0 + 1, ny - 1);

    const varType f00 = Get(i0, j0);
    const varType f10 = Get(i1, j0);
    const varType f01 = Get(i0, j1);
    const varType f11 = Get(i1, j1);

    assert(std::isfinite((REAL_LITERAL(1.0) - fy) *
                             ((REAL_LITERAL(1.0) - fx) * f00 + fx * f10) +
                         fy * ((REAL_LITERAL(1.0) - fx) * f01 + fx * f11)));
    return (REAL_LITERAL(1.0) - fy) *
               ((REAL_LITERAL(1.0) - fx) * f00 + fx * f10) +
           fy * ((REAL_LITERAL(1.0) - fx) * f01 + fx * f11);
  }

  varType *A = nullptr;
};
