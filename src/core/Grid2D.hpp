#pragma once
#include "Precision.hpp"
#include <string>
#include <vector>

class Grid2D {
public:
  size_t nx, ny;
  varType dx, dy;
  std::vector<varType> A;

  Grid2D(size_t nx, size_t ny) : nx(nx), ny(ny), A(nx * ny, 0.0) {}

  // manipulating grid values
  varType Get(size_t i, size_t j) const;
  void Set(size_t i, size_t j, varType val);

  // utility functions
  bool InBounds(size_t i, size_t j);

  varType Interpolate(varType x, varType y,varType dx, varType dy);
  
  void InitRectangle(varType constVel);
};
