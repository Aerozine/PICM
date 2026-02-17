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
      // int xCenter = std::abs(midX - i);
      this->Set(i, j, constVel);
    }
  }
  return;
}

varType Grid2D::Interpolate(varType x, varType y, varType dx, varType dy,
                            int field) {
  // if (x <= 0 || x >= this->nx * dx) assert(false);
  // if (y <= 0 || y >= this->ny * dy) assert(false);

  // if (field == 0) y -= dy / 2; // handle field u
  // if (field == 1) x -= dx / 2; // handle field v

  varType k = x / dx;
  varType l = y / dy;

  if (field == 0) { // u : décalé de +0.5 en y
    l -= varType(0.5);
  } else if (field == 1) { // v : décalé de +0.5 en x
    k -= varType(0.5);
  } else {
  }

  int i0 = std::floor(k);
  int j0 = std::floor(l);

  // printf("i0 = %d\n", i0);
  // printf("j0 = %d\n", j0);

  varType a = k - (varType)i0;
  varType b = l - (varType)j0;

  // printf("a = %f\n", a);
  // printf("b = %f\n", b);

  varType f00, f10, f01, f11;

  if (k < 0) { // v does not have 4 corners to interpolate from
    // printf("entered k case \n");
    f00 = 0.0;
    f01 = 0.0; // impose homogeneous BC on ghost cells

    varType i1 = std::ceil(k);
    f10 = Get(i1, j0);
    f11 = Get(i1, j0 + 1);
  }

  else if (l < 0) { // u does not have 4 corners to interpolate from

    // printf("entered l case \n");
    f00 = 0.0;
    f10 = 0.0; // impose homogeneous BC on ghost cells

    varType j1 = std::ceil(l);
    f01 = Get(i0, j1);
    f11 = Get(i0 + 1, j1);
  }

  else {
    // printf("entered default case \n");
    f00 = Get(i0, j0);
    // printf("f00 = %f\n", f00);
    f10 = Get(i0 + 1, j0);
    // printf("f10 = %f\n", f10);
    f01 = Get(i0, j0 + 1);
    // printf("f01 = %f\n", f01);
    f11 = Get(i0 + 1, j0 + 1);
    // printf("f11 = %f\n", f11);
  }

  return (1 - a) * (1 - b) * f00 + a * (1 - b) * f10 + (1 - a) * b * f01 +
         a * b * f11;
}
