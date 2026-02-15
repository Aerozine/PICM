#include "Grid2D.hpp"
#include <cassert>
#include <random>

bool Grid2D::InBounds(int i, int j) const { return (i < nx) && (j < ny); }

void Grid2D::InitRectangle(const varType constVel) {
  const int midX = nx / 2;
  const int midY = ny / 2;

  constexpr int offsetX = 10;
  const int offsetY = 10;

  for (int i = midX - offsetX; i < midX + offsetX; i++) {
    for (int j = midY - offsetY; j < midY + offsetY; j++) {
      //int xCenter = std::abs(midX - i);
      this->Set(i, j, constVel);
    }
  }
}

varType Grid2D::Interpolate(varType x, varType y,varType dx, varType dy) const
{
    if (x <= 0 || x >= (this->nx - 1) * dx) assert(false);
    if (y <= 0 || y >= (this->ny - 1) * dy) assert(false);

    const varType k = x / dx;
    const varType l = y / dy;

    const int i0 = std::floor(k);
    const int j0 = std::floor(l);

    const varType a = k - static_cast<varType>(i0);
    const varType b = l - static_cast<varType>(j0);

    const int i = (int) i0;
    const int j = (int) j0;

    const varType f00 = Get(i, j);
    const varType f10 = Get(i + 1, j);
    const varType f01 = Get(i, j + 1);
    const varType f11 = Get(i + 1, j + 1);

    return (1 - a) * (1 - b) * f00
          + a * (1 - b) * f10
          + (1 - a) * b * f01
          + a * b * f11;
}
