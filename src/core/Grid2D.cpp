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

Grid2D::varType Grid2D::InterpolateCenter(size_t iCenter, size_t jCenter,
                                          const size_t field) {
  // iCenter/jCenter = numbering of pressure grid
  varType interpol;
 
  switch (field) {
    case 0: interpol = (this -> Get(iCenter, jCenter) 
                        + this -> Get(iCenter, jCenter + 1)) / 2;   
              break;

    case 1: interpol = (this -> Get(iCenter, jCenter)   
                        + this -> Get(iCenter + 1, jCenter)) / 2;   
              break;

    default: assert(false); 
             break;
  }
  return interpol;
}


