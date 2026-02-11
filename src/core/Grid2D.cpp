#include "Grid2D.hpp"
#include <cassert>
#include <random>
#include <stdio.h>

varType Grid2D::Get(size_t i, size_t j) const { return A[nx * j + i]; }

void Grid2D::Set(size_t i, size_t j, varType val) {
  A[nx * j + i] = val;
  return;
}

bool Grid2D::InBounds(size_t i, size_t j) { return (i < nx) && (j < ny); }

void Grid2D::InitRectangle(varType constVel) {
  int midX = nx / 2;
  int midY = ny / 2;

  int offsetX = 10;
  int offsetY = 10;

  for (int i = midX - offsetX; i < midX + offsetX; i++) {
    for (int j = midY - offsetY; j < midY + offsetY; j++) {
      //int xCenter = std::abs(midX - i);
      this->Set(i, j, constVel);
    }
  }
  return;
}

varType Grid2D::Interpolate(varType x, varType y,varType dx, varType dy) {
    if (x <= 0 || x >= (this->nx - 1) * dx) assert(false);
    if (y <= 0 || y >= (this->ny - 1) * dy) assert(false);

    varType k = x / dx;
    varType l = y / dy;

    int i0 = std::floor(k);
    int j0 = std::floor(l);

    varType a = k - (varType) i0; 
    varType b = l - (varType) j0;  

    int i = (size_t) i0;
    int j = (size_t) j0;

    varType f00 = Get(i, j);
    varType f10 = Get(i + 1, j);
    varType f01 = Get(i, j + 1);
    varType f11 = Get(i + 1, j + 1);

    return (1 - a) * (1 - b) * f00
          + a * (1 - b) * f10
          + (1 - a) * b * f01
          + a * b * f11;
}
