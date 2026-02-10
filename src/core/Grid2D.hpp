#pragma once
#include "Precision.hpp"
#include <string>
#include <vector>

class Grid2D {
public:
  size_t nx, ny;
  Real dx, dy;
  std::vector<Real> A;

  Grid2D(size_t nx, size_t ny) : nx(nx), ny(ny), A(nx * ny, 0.0) {}

  // manipulating grid values
  Real Get(size_t i, size_t j) const { return A[nx * j + i]; }
  void Set(size_t i, size_t j, Real val) {A[nx * j + i] = val;}

  // utility functions
  bool InBounds(size_t i, size_t j);
  void FillRandom();
  Real InterpolateCenter(size_t iCenter, size_t jCenter, const size_t field);
};

