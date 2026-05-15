#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>

#include "Precision.hpp"

///@brief  **ROW-MAJOR 2d grid** @c A[nx*j+i]
class Grid2D {
public:
  int nx = 0;
  int ny = 0;

  Grid2D() = default;

  Grid2D(int nx_, int ny_)
      : nx(nx_), ny(ny_),
        A(new varType[static_cast<std::size_t>(nx_) * ny_]()) {}

  ~Grid2D() { delete[] A; }

  Grid2D(const Grid2D &o)
      : nx(o.nx), ny(o.ny),
        A(new varType[static_cast<std::size_t>(o.nx) * o.ny]) {
    std::copy(o.A, o.A + static_cast<std::size_t>(nx) * ny, A);
  }

  Grid2D &operator=(const Grid2D &o) {
    if (this == &o)
      return *this;
    if (nx * ny != o.nx * o.ny) {
      delete[] A;
      A = new varType[static_cast<std::size_t>(o.nx) * o.ny];
    }
    nx = o.nx;
    ny = o.ny;
    std::copy(o.A, o.A + static_cast<std::size_t>(nx) * ny, A);
    return *this;
  }

  Grid2D(Grid2D &&o) noexcept : nx(o.nx), ny(o.ny), A(o.A) {
    o.nx = 0;
    o.ny = 0;
    o.A = nullptr;
  }

  Grid2D &operator=(Grid2D &&o) noexcept {
    if (this == &o)
      return *this;
    delete[] A;
    nx = o.nx;
    ny = o.ny;
    A = o.A;
    o.nx = 0;
    o.ny = 0;
    o.A = nullptr;
    return *this;
  }

  void reset() noexcept {
    std::fill(A, A + static_cast<std::size_t>(nx) * ny, varType{});
  }

  [[nodiscard]] inline varType &operator()(int i, int j) noexcept {
    return A[nx * j + i];
  }
  [[nodiscard]] inline varType operator()(int i, int j) const noexcept {
    return A[nx * j + i];
  }

  /// @brief Read the value at cell (i, j).
  [[nodiscard]] varType Get(const int i, const int j) const {
    return A[nx * j + i];
  }

  /// @brief Write a value into cell (i, j).
  void Set(const int i, const int j, const varType val) { A[nx * j + i] = val; }

  varType *A = nullptr;
};
