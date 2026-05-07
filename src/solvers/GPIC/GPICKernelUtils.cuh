#pragma once
// Device-only utilities for GPIC kernels.
// Included only by GPIC.cu (compiled by nvcc).

#include "../../core/Fields.hpp"
#include "../../core/Precision.hpp"

// ─── Kernel weight (mirrors PIC::hat) ────────────────────────────────────────

__device__ __forceinline__ varType hat_device(varType r) {
#ifdef USE_SPEED
    const varType a = (r < varType(0)) ? -r : r;
    return (a < varType(1)) ? (varType(1) - a) : varType(0);
#else
    if (r >= varType(-1.5) && r < varType(-0.5))
        return varType(0.5) * (r + varType(1.5)) * (r + varType(1.5));
    if (r >= varType(-0.5) && r < varType(0.5))
        return varType(0.75) - r * r;
    if (r >= varType(0.5) && r < varType(1.5))
        return varType(0.5) * (varType(1.5) - r) * (varType(1.5) - r);
    return varType(0);
#endif
}

// ─── Bilinear interpolation (mirrors Grid2D::interpolate<field>) ──────────────
// gnx / gny are the Grid2D::nx and Grid2D::ny members.
// field == 0 → u-face (staggered in y), field == 1 → v-face (staggered in x).

template <unsigned char field>
__device__ __forceinline__ varType bilinear_device(
    varType x, varType y,
    varType dx, varType dy,
    const varType* __restrict__ grid,
    int gnx, int gny)
{
    varType i_real = x / dx;
    varType j_real = y / dy;

    if (field == 0) j_real -= varType(0.5);
    if (field == 1) i_real -= varType(0.5);

    if (i_real < varType(0))       i_real = varType(0);
    if (i_real > varType(gnx - 1)) i_real = varType(gnx - 1);
    if (j_real < varType(0))       j_real = varType(0);
    if (j_real > varType(gny - 1)) j_real = varType(gny - 1);

    const int i0 = (int)i_real;
    const int j0 = (int)j_real;
    const int i1 = (i0 + 1 < gnx) ? i0 + 1 : gnx - 1;
    const int j1 = (j0 + 1 < gny) ? j0 + 1 : gny - 1;

    const varType fx = i_real - varType(i0);
    const varType fy = j_real - varType(j0);

    const varType f00 = grid[gnx * j0 + i0];
    const varType f10 = grid[gnx * j0 + i1];
    const varType f01 = grid[gnx * j1 + i0];
    const varType f11 = grid[gnx * j1 + i1];

    return (varType(1) - fy) * ((varType(1) - fx) * f00 + fx * f10)
         + fy               * ((varType(1) - fx) * f01 + fx * f11);
}

// ─── Domain / solid checks ───────────────────────────────────────────────────

// Returns true if (x,y) is inside the domain and not in a solid cell.
// d_labels is the (nx+2)*(ny+2) label array; physical cell (ci,cj) is at
// label index (nx+2)*(cj+1) + (ci+1).
__device__ __forceinline__ bool isFluidPoint_device(
    varType x, varType y,
    int nx, int ny, varType dx, varType dy,
    const uint16_t* __restrict__ d_labels)
{
    if (x < varType(0) || x >= dx * nx || y < varType(0) || y >= dy * ny)
        return false;
    int ci = (int)(x / dx);
    int cj = (int)(y / dy);
    if (ci >= nx) ci = nx - 1;
    if (cj >= ny) cj = ny - 1;
    const uint16_t lbl = d_labels[(nx + 2) * (cj + 1) + (ci + 1)];
    return (lbl & Fields2D::SOLID) == 0;
}

// Binary bisection (12 iterations) to find the last non-solid point on the
// segment from (x0,y0) to (x1,y1). Mirrors PIC::Advect clampToLastFluidPoint.
__device__ __forceinline__ void clampToLastFluidPoint_device(
    varType x0, varType y0, varType x1, varType y1,
    varType& xc, varType& yc,
    int nx, int ny, varType dx, varType dy,
    const uint16_t* __restrict__ d_labels)
{
    varType lo = varType(0);
    varType hi = varType(1);
    for (int it = 0; it < 12; ++it) {
        const varType mid = varType(0.5) * (lo + hi);
        const varType xm  = x0 + mid * (x1 - x0);
        const varType ym  = y0 + mid * (y1 - y0);
        if (isFluidPoint_device(xm, ym, nx, ny, dx, dy, d_labels))
            lo = mid;
        else
            hi = mid;
    }
    const varType backed = (lo > REAL_LITERAL(1e-4)) ? lo - REAL_LITERAL(1e-4) : varType(0);
    xc = x0 + backed * (x1 - x0);
    yc = y0 + backed * (y1 - y0);
    if (!isFluidPoint_device(xc, yc, nx, ny, dx, dy, d_labels)) {
        xc = x0;
        yc = y0;
    }
}
