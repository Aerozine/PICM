#include "Grid2D.hpp"
#include <cassert>
#include <random>
#include <stdio.h>

varType Grid2D::Get(int i, int j) const { return A[ny * i + j]; }

void Grid2D::Set(int i, int j, varType val) {
  A[ny * i + j] = val;
  return;
}

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
    if (x <= 0 || x >= this->nx * dx) assert(false);
    if (y <= 0 || y >= this->ny * dy) assert(false);

    varType k = x / dx;
    varType l = y / dy;

    int i0 = std::floor(k);
    int j0 = std::floor(l);

    varType a = k - (varType) i0; 
    varType b = l - (varType) j0;  

    varType f00 = Get(i0, j0);
    varType f10 = Get(i0 + 1, j0);
    varType f01 = Get(i0, j0 + 1);
    varType f11 = Get(i0 + 1, j0 + 1);

    return (1 - a) * (1 - b) * f00
          + a * (1 - b) * f10
          + (1 - a) * b * f01
          + a * b * f11;
}

