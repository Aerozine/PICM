#include "Grid2D.hpp"
#include <random>
#include <cassert>

// prefer some macros ?
Grid2D::varType Grid2D::Get(size_t i, size_t j) const{
  return A[nx * j + i];
}

void Grid2D::Set(size_t i, size_t j, varType val){
  A[nx * j + i] = val;
  return;
}

// utility functions
bool Grid2D::InBounds(size_t i, size_t j){
  return (i < nx) && (j < ny);
}

void Grid2D::InitRectangle(varType constVel = 10) {
  int midX = nx/2;
  int midY = ny/2;

  int offsetX = 20;
  int offsetY = 10;

  for(int i = midX - offsetX; i < midX + offsetX; i++){
    for(int j = midY - offsetY; j < midY + offsetY; j++){
      this -> Set(i, j, constVel); 
    }
  }
  return;
}


void Grid2D::FillRandom() {
  std::mt19937 gen(42); // fixed seed
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  for(size_t j = 0; j < ny; j++){
    for(size_t i = 0; i < nx; i++){
      this->Set(i, j, dist(gen));
    }
  }
  return;
}

/*
Grid2D::varType Grid2D::Interpolate(varType x, varType y, varType dx, varType dy){
  size_t i = std::floor(x / dx);
  size_t j = std::floor(y / dy);

  varType restX = x % i;
  varType restY = y % j;
  
  return 0.0;
}
*/

Grid2D::varType Grid2D::Interpolate(varType sx , varType sy,
                                    varType dx, varType dy, size_t field) {

  if (sx >= nx && sy >= ny) assert(false); 
  
  // Determine the region [xi, xi+1]
  size_t si = std::floor(sx / dx);
  size_t sj = std::floor(sy / dx);

  varType interpol; 
 
  switch (field) {
    case 0: { 
              varType alpha = (sx - si) / dx;
              interpol = (1 - alpha) * (this -> Get(si, sj)) 
                     + alpha * (this -> Get(si, sj + 1));
            }
            break;
    
    case 1: { 
              varType alpha = (sy - si) / dy;
              interpol = (1 - alpha) * (this -> Get(si, sj))
                     + alpha * (this -> Get(si + 1, sj));
            }
            break;
    
    default: assert(false); 
             break;
  }
  return interpol;
}

