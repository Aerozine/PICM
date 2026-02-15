#pragma once
#include "Precision.hpp"
#include <string>
#include <vector>

class Grid2D {
public:
  int nx, ny;
  //varType dx, dy;
  std::vector<varType> A;

  Grid2D(const int nx, const int ny) : nx(nx), ny(ny), A(nx * ny, 0.0) {}

  // Inline functions defined directly in the class
  [[nodiscard]] varType Get(int i, int j) const {
    return A[nx * j + i]; 
  }

  void Set(int i, int j, const varType val) {
    A[nx * j + i] = val;
  }

  // utility functions
  [[nodiscard]] bool InBounds(int i, int j) const;

  [[nodiscard]] varType Interpolate(varType x, varType y,varType dx, varType dy) const;
  
  void InitRectangle(varType constVel);
};
