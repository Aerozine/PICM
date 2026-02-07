#include "Grid2D.hpp"

// prefer some macros ?
Grid2D::varType Grid2D::Get(size_t i, size_t j) const{
  return A[nx * (j) + (i)];
}

void Grid2D::Set(size_t i, size_t j, varType val){
  A[nx * (j) + (i)] = val;
  return;
}

// utility functions
Grid2D Grid2D::InitRandomGrid(size_t nx, size_t ny) {
  Grid2D grid(nx, ny);

  std::mt19937 gen(42); // fixed seed
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  for(size_t j = 0; j < ny; j++){
    for(size_t i = 0; i < nx; i++){
      grid.Set(i, j, dist(gen));
    }
  }
  return grid;
}

Grid2D::varType Grid2D::Interpolate(size_t nx, size_t ny){
  return 0.0;
}
