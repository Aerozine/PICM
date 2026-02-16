#pragma once
#include "Precision.hpp"
#include <string>
#include <vector>

class Grid2D {
public:
  int nx, ny;
  varType dx, dy;
  std::vector<varType> A;

  Grid2D(int nx, int ny) : nx(nx), ny(ny), A(nx * ny, 0.0) {}

  // manipulating grid values
  varType Get(int i, int j) const;
  void Set(int i, int j, varType val);

  // utility functions
  bool InBounds(int i, int j) { return (i < nx) && (j < ny); }
  void InitRectangle(varType constVel);
  varType Interpolate(varType x, varType y, varType dx, varType dy, int field);
};
