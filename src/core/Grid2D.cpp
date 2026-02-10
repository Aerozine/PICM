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

  int offsetX = 4;
  int offsetY = 4;

  for (int i = midX - offsetX; i < midX + offsetX; i++) {
    for (int j = midY - offsetY; j < midY + offsetY; j++) {
      //int xCenter = std::abs(midX - i);
      this->Set(i, j, constVel);
    }
  }
  return;
}

// needs to be redone with bilinear interpolation
varType Grid2D::Interpolate(varType sx, varType sy, varType dx, varType dy,
                            size_t field) {

  if (sx >= nx && sy >= ny)
    assert(false);

  // Determine the region [xi, xi+1]
  size_t si = std::floor(sx / dx);
  size_t sj = std::floor(sy / dx);

  varType interpol;

  switch (field) {
  case 0: {
    varType alpha = (sx - si) / dx;
    interpol =
        (1 - alpha) * (this->Get(si, sj)) + alpha * (this->Get(si, sj + 1));
  } break;

  case 1: {
    varType alpha = (sy - si) / dy;
    interpol =
        (1 - alpha) * (this->Get(si, sj)) + alpha * (this->Get(si + 1, sj));
  } break;

  default:
    assert(false);
    break;
  }
  return interpol;
}
