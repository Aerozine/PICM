#pragma once
#include <vector>

typedef float varType;

class Grid2D {
public:
  size_t nx, ny;
  std::vector<std::vector<varType>> A;
  
  Grid2D(size_t nx, size_t ny)
    : nx(nx), ny(ny){}
  
  void InitRandomGrid(size_t nx, size_t ny); 
  Grid2D Div(Grid2D grid);
  varType Interpolate(Grid2D grid, size_t nx, size_t ny);
};
