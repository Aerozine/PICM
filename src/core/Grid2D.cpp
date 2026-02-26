#include "Grid2D.hpp"
#include <algorithm>
#include <cmath>

// Bilinear interpolation
//
// Staggered MAC grid offsets:
//
//   field 0 (u): nodes at (i·dx,       (j+0.5)·dy)  →  j_real -= 0.5
//   field 1 (v): nodes at ((i+0.5)·dx,  j·dy      )  →  i_real -= 0.5
//   other (p, cell-centred): no offset.
//
// After applying the offset, (i_real, j_real) is the continuous index into
// the node array. We split it into an integer base (i0, j0) and a fractional
// weight (fx, fy), clamp the base so the 2×2 stencil stays in bounds, then
// bilinearly blend the four surrounding values.
//
// Get()/Set() use row-major storage (A[nx*j + i]) so the index arithmetic
// is encapsulated there — this function is unchanged relative to the
// column-major version.

varType Grid2D::Interpolate(varType x, varType y, varType dx, varType dy,
                            int field) const {
  varType i_real = x / dx;
  varType j_real = y / dy;

  if (field == 0)
    j_real -= REAL_LITERAL(0.5); // u-face: staggered in y
  else if (field == 1)
    i_real -= REAL_LITERAL(0.5); // v-face: staggered in x

  int i0 = static_cast<int>(std::floor(i_real));
  int j0 = static_cast<int>(std::floor(j_real));

  const varType fx = i_real - static_cast<varType>(i0);
  const varType fy = j_real - static_cast<varType>(j0);

  // Clamp so that (i0, j0) … (i0+1, j0+1) are all valid indices.
  i0 = std::clamp(i0, 0, nx - 2);
  j0 = std::clamp(j0, 0, ny - 2);

  const varType f00 = Get(i0,     j0    );
  const varType f10 = Get(i0 + 1, j0    );
  const varType f01 = Get(i0,     j0 + 1);
  const varType f11 = Get(i0 + 1, j0 + 1);

  return (REAL_LITERAL(1.0) - fy) *
             ((REAL_LITERAL(1.0) - fx) * f00 + fx * f10) +
         fy * ((REAL_LITERAL(1.0) - fx) * f01 + fx * f11);
}
