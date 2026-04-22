#include "Grid2D.hpp"
#include <algorithm>
#include <cmath>

varType Grid2D::interpolate(const varType x, const varType y, const varType dx,
                            const varType dy, const int field) const {
  varType i_real = x / dx;
  varType j_real = y / dy;

  if (field == 0)
    j_real -= REAL_LITERAL(0.5); // u-face: staggered in y
  else if (field == 1)
    i_real -= REAL_LITERAL(0.5); // v-face: staggered in x
  // else no staggering

  int i0 = static_cast<int>(std::floor(i_real));
  int j0 = static_cast<int>(std::floor(j_real));

  const varType fx = i_real - static_cast<varType>(i0);
  const varType fy = j_real - static_cast<varType>(j0);

  i0 = std::clamp(i0, 0, nx - 1);
  j0 = std::clamp(j0, 0, ny - 1);

  const varType f00 = Get(i0, j0);
  const varType f10 = Get(i0 + 1, j0);
  const varType f01 = Get(i0, j0 + 1);
  const varType f11 = Get(i0 + 1, j0 + 1);

  assert(std::isfinite(
  (REAL_LITERAL(1.0) - fy) *
           ((REAL_LITERAL(1.0) - fx) * f00 + fx * f10) +
       fy * ((REAL_LITERAL(1.0) - fx) * f01 + fx * f11)
       ));
  return (REAL_LITERAL(1.0) - fy) *
             ((REAL_LITERAL(1.0) - fx) * f00 + fx * f10) +
         fy * ((REAL_LITERAL(1.0) - fx) * f01 + fx * f11);
}
