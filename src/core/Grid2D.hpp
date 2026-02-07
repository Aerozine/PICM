#pragma once
#include <vector>

class Grid2D {
public:

  typedef float varType;

  size_t nx, ny;
  std::vector<varType> A;
  
  Grid2D(size_t nx, size_t ny)
    : nx(nx), ny(ny),
      A(nx * ny, 0.0) {}
  
  // manipulating grid values 
  varType Get(size_t i, size_t j) const;
  void Set(size_t i, size_t j, varType val);

  // utility functions
  static Grid2D InitRandomGrid(size_t nx, size_t ny); 
  Grid2D Div();
  varType Interpolate(size_t nx, size_t ny);
};
