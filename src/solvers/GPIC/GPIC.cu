#ifdef USE_CUDA

#include "GPIC.hpp"
#include "GPICKernelUtils.cuh"

#include <cuda/functional>
#include <cuda_runtime.h>
#include <cuda/std/functional>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

// ─── Error checking ───────────────────────────────────────────────────────────

#define CUDA_CHECK(call)                                                      \
  do {                                                                        \
    cudaError_t _e = (call);                                                  \
    if (_e != cudaSuccess) {                                                  \
      std::fprintf(stderr, "[CUDA] %s:%d  %s\n", __FILE__, __LINE__,         \
                   cudaGetErrorString(_e));                                   \
      std::exit(EXIT_FAILURE);                                                \
    }                                                                         \
  } while (0)

// ─── DeviceState ─────────────────────────────────────────────────────────────

struct GPIC::DeviceState {
    // Flat SoA particle arrays  [maxParticles]
    varType* d_pos_x    = nullptr;
    varType* d_pos_y    = nullptr;
    varType* d_vel_x    = nullptr;
    varType* d_vel_y    = nullptr;
    varType* d_cu_x     = nullptr;
    varType* d_cu_y     = nullptr;
    varType* d_cv_x     = nullptr;
    varType* d_cv_y     = nullptr;
    int*     d_cell_idx = nullptr;
    int*     d_alive    = nullptr;

    // Grid velocity fields (persistent on device between substeps)
    varType*  d_u      = nullptr;  // (nx+1) * ny
    varType*  d_v      = nullptr;  // nx * (ny+1)
    varType*  d_u_old  = nullptr;  // FLIP transfer cache
    varType*  d_v_old  = nullptr;

    // Pressure + CG workspace (all persistent = warm start between steps)
    varType*  d_p         = nullptr;  // (nx+2) * (ny+2)  pressure (warm-start across steps)
    varType*  d_rhs       = nullptr;  // (nx+2) * (ny+2)  RHS  b = -coef * div(u,v)
    varType*  d_cg_r      = nullptr;  // (nx+2) * (ny+2)  CG residual  r = b - A*p
    varType*  d_cg_d      = nullptr;  // (nx+2) * (ny+2)  CG search direction
    varType*  d_cg_q      = nullptr;  // (nx+2) * (ny+2)  CG matrix-vector product  A*d
    // CG scalar state kept on device to avoid per-iteration D2H syncs
    // Layout: [0]=sigma (r·r), [1]=d·Ad, [2]=sigma_new (r_new·r_new)
    double*   d_cg_scalars = nullptr;  // 3 doubles

    // P2G accumulators (cleared each substep)
    varType*  d_u_num  = nullptr;  // (nx+1) * ny
    varType*  d_u_den  = nullptr;
    varType*  d_v_num  = nullptr;  // nx * (ny+1)
    varType*  d_v_den  = nullptr;

    // Labels (persistent; updated by cellLabelKernel + airCleanupKernel)
    uint16_t* d_labels = nullptr;  // (nx+2) * (ny+2)

    // Per-cell particle count (rebuilt each step via countCellsKernel)
    int* d_cell_count  = nullptr;  // [nx*ny]

    // Reusable scratch (avoid per-call alloc)
    varType* d_tmp     = nullptr;  // [maxParticles] varType scratch
    int*     d_perm    = nullptr;  // [maxParticles] sort permutation
    int*     d_particle_count = nullptr;   // device-side refill counter
    int*     d_refill_overflow = nullptr;  // skipped particles when maxParticles is exceeded

    int maxParticles = 0;
    int nParticles   = 0;
    int nx = 0, ny = 0;
    varType dx{}, dy{}, dt{};
    varType density{};
    int refillSeed = 0;
    int transferMode = 0; // 0=PIC, 1=FLIP, 2=APIC
    bool needsAffine = false;
};

// ─── Functors ─────────────────────────────────────────────────────────────────

struct NonZero {
    __host__ __device__ bool operator()(int x) const { return x != 0; }
};
struct SquareFunctor {
    __host__ __device__ varType operator()(varType x) const { return x * x; }
};
struct AbsFunctor {
    __host__ __device__ varType operator()(varType x) const {
        return x < varType(0) ? -x : x;
    }
};

static constexpr int BLOCK     = 256;
static constexpr int BX        = 32;
static constexpr int BY        = 8;
static constexpr int DOT_GRIDS = 64;   // blocks for dot-product reductions
static constexpr int CG_CHECK  = 8;    // D2H convergence check every N CG iters
static constexpr int TRANSFER_PIC  = 0;
static constexpr int TRANSFER_FLIP = 1;
static constexpr int TRANSFER_APIC = 2;

// ─── P2G scatter ──────────────────────────────────────────────────────────────

__global__ void p2gKernel(
    int nParticles, int nx, int ny, int transferMode,
    varType dx, varType dy,
    const varType* __restrict__ d_pos_x,
    const varType* __restrict__ d_pos_y,
    const varType* __restrict__ d_vel_x,
    const varType* __restrict__ d_vel_y,
    const varType* __restrict__ d_cu_x,
    const varType* __restrict__ d_cu_y,
    const varType* __restrict__ d_cv_x,
    const varType* __restrict__ d_cv_y,
    varType* d_u_num, varType* d_u_den,
    varType* d_v_num, varType* d_v_den)
{
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= nParticles) return;

    const varType xg  = d_pos_x[p] / dx;
    const varType yg  = d_pos_y[p] / dy;
    const varType vx  = d_vel_x[p];
    const varType vy  = d_vel_y[p];

    const int ci_p = min((int)xg, nx - 1);
    const int cj_p = min((int)yg, ny - 1);

#ifdef USE_SPEED
    const int u_fi_lo = max(0,    ci_p - 1);
    const int u_fi_hi = min(nx,   ci_p + 2);
    const int u_fj_lo = max(0,    cj_p - 2);
    const int u_fj_hi = min(ny-1, cj_p + 1);
#else
    const int u_fi_lo = max(0,    ci_p - 2);
    const int u_fi_hi = min(nx,   ci_p + 3);
    const int u_fj_lo = max(0,    cj_p - 3);
    const int u_fj_hi = min(ny-1, cj_p + 2);
#endif
    for (int fj = u_fj_lo; fj <= u_fj_hi; ++fj) {
        for (int fi = u_fi_lo; fi <= u_fi_hi; ++fi) {
            const varType w = hat_device(xg - varType(fi))
                            * hat_device(yg - varType(fj) - varType(0.5));
            if (w > varType(0)) {
                const int id = (nx + 1) * fj + fi;
                const varType faceX = varType(fi) * dx;
                const varType faceY = (varType(fj) + varType(0.5)) * dy;
                const varType value =
                    (transferMode == TRANSFER_APIC)
                        ? vx + d_cu_x[p] * (faceX - d_pos_x[p])
                             + d_cu_y[p] * (faceY - d_pos_y[p])
                        : vx;
                atomicAdd(&d_u_num[id], w * value);
                atomicAdd(&d_u_den[id], w);
            }
        }
    }

#ifdef USE_SPEED
    const int v_fi_lo = max(0,    ci_p - 2);
    const int v_fi_hi = min(nx-1, ci_p + 1);
    const int v_fj_lo = max(0,    cj_p - 1);
    const int v_fj_hi = min(ny,   cj_p + 2);
#else
    const int v_fi_lo = max(0,    ci_p - 3);
    const int v_fi_hi = min(nx-1, ci_p + 2);
    const int v_fj_lo = max(0,    cj_p - 2);
    const int v_fj_hi = min(ny,   cj_p + 3);
#endif
    for (int fj = v_fj_lo; fj <= v_fj_hi; ++fj) {
        for (int fi = v_fi_lo; fi <= v_fi_hi; ++fi) {
            const varType w = hat_device(xg - varType(fi) - varType(0.5))
                            * hat_device(yg - varType(fj));
            if (w > varType(0)) {
                const int id = nx * fj + fi;
                const varType faceX = (varType(fi) + varType(0.5)) * dx;
                const varType faceY = varType(fj) * dy;
                const varType value =
                    (transferMode == TRANSFER_APIC)
                        ? vy + d_cv_x[p] * (faceX - d_pos_x[p])
                             + d_cv_y[p] * (faceY - d_pos_y[p])
                        : vy;
                atomicAdd(&d_v_num[id], w * value);
                atomicAdd(&d_v_den[id], w);
            }
        }
    }
}

// ─── P2G normalise ────────────────────────────────────────────────────────────

__global__ void p2gNormaliseKernel_U(
    int nx, int ny,
    const varType* __restrict__ d_u_num,
    const varType* __restrict__ d_u_den,
    const uint16_t* __restrict__ d_labels,
    varType* d_u)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= (nx + 1) * ny) return;

    const int fi = id % (nx + 1);
    const int fj = id / (nx + 1);

    const uint16_t lbl_L = d_labels[(nx + 2) * (fj + 1) + fi];
    const uint16_t lbl_R = d_labels[(nx + 2) * (fj + 1) + fi + 1];

    if (lbl_L & Fields2D::BC_U) return;

    if ((lbl_L & Fields2D::SOLID) || (lbl_R & Fields2D::SOLID)) {
        d_u[id] = varType(0);
        return;
    }

    const varType den = d_u_den[id];
    d_u[id] = (den >= REAL_EPSILON) ? d_u_num[id] / den : varType(0);
}

__global__ void p2gNormaliseKernel_V(
    int nx, int ny,
    const varType* __restrict__ d_v_num,
    const varType* __restrict__ d_v_den,
    const uint16_t* __restrict__ d_labels,
    varType* d_v)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= nx * (ny + 1)) return;

    const int fi = id % nx;
    const int fj = id / nx;

    const uint16_t lbl_B = d_labels[(nx + 2) * fj       + fi + 1];
    const uint16_t lbl_T = d_labels[(nx + 2) * (fj + 1) + fi + 1];

    if (lbl_B & Fields2D::BC_V) return;

    if ((lbl_B & Fields2D::SOLID) || (lbl_T & Fields2D::SOLID)) {
        d_v[id] = varType(0);
        return;
    }

    const varType den = d_v_den[id];
    d_v[id] = (den >= REAL_EPSILON) ? d_v_num[id] / den : varType(0);
}

// ─── G2P ─────────────────────────────────────────────────────────────────────

__device__ __forceinline__ void accumulateAffineDevice(
    const varType* __restrict__ grid,
    varType xg, varType yg,
    int gnx, int gny,
    varType dx, varType dy,
    int radius,
    varType& value,
    varType& gradX,
    varType& gradY)
{
    const int i0 = static_cast<int>(floor(xg));
    const int j0 = static_cast<int>(floor(yg));

    varType valueRaw = varType(0);
    varType gradXRaw = varType(0);
    varType gradYRaw = varType(0);
    varType weightSum = varType(0);
    varType dWeightX = varType(0);
    varType dWeightY = varType(0);

    for (int dj = -radius; dj <= radius; ++dj) {
        for (int di = -radius; di <= radius; ++di) {
            const int i = i0 + di;
            const int j = j0 + dj;
            if (i < 0 || i >= gnx || j < 0 || j >= gny) continue;

            const varType wx = hat_device(xg - varType(i));
            const varType wy = hat_device(yg - varType(j));
            const varType w = wx * wy;
            if (w <= varType(0)) continue;

            const varType sample = grid[gnx * j + i];
            const varType dwx = (dhat_device(xg - varType(i)) / dx) * wy;
            const varType dwy = wx * (dhat_device(yg - varType(j)) / dy);

            valueRaw += w * sample;
            gradXRaw += dwx * sample;
            gradYRaw += dwy * sample;
            weightSum += w;
            dWeightX += dwx;
            dWeightY += dwy;
        }
    }

    if (weightSum <= REAL_EPSILON) {
        value = varType(0);
        gradX = varType(0);
        gradY = varType(0);
        return;
    }

    value = valueRaw / weightSum;
    const varType invWeightSq = varType(1) / (weightSum * weightSum);
    gradX = (gradXRaw * weightSum - valueRaw * dWeightX) * invWeightSq;
    gradY = (gradYRaw * weightSum - valueRaw * dWeightY) * invWeightSq;
}

__global__ void g2pKernel(
    int nParticles, int nx, int ny, int transferMode, int kernelRadius,
    varType dx, varType dy, varType dt, varType gravity, varType coefPic,
    const varType* __restrict__ d_u,
    const varType* __restrict__ d_v,
    const varType* __restrict__ d_u_old,
    const varType* __restrict__ d_v_old,
    const varType* __restrict__ d_pos_x,
    const varType* __restrict__ d_pos_y,
    varType* d_vel_x,
    varType* d_vel_y,
    varType* d_cu_x,
    varType* d_cu_y,
    varType* d_cv_x,
    varType* d_cv_y)
{
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= nParticles) return;

    const varType x = d_pos_x[p];
    const varType y = d_pos_y[p];

    const varType uNew = bilinear_device<0>(x, y, dx, dy, d_u, nx + 1, ny);
    const varType vNew = bilinear_device<1>(x, y, dx, dy, d_v, nx, ny + 1);

    if (transferMode == TRANSFER_FLIP) {
        const varType coefFlip = varType(1) - coefPic;
        const varType uOld = bilinear_device<0>(x, y, dx, dy, d_u_old, nx + 1, ny);
        const varType vOld = bilinear_device<1>(x, y, dx, dy, d_v_old, nx, ny + 1);
        d_vel_x[p] = coefPic * uNew + coefFlip * (d_vel_x[p] + (uNew - uOld));
        d_vel_y[p] = coefPic * vNew + coefFlip * (d_vel_y[p] + (vNew - vOld))
                   - dt * gravity;
        return;
    }

    if (transferMode == TRANSFER_APIC) {
        varType u = varType(0), cuX = varType(0), cuY = varType(0);
        varType v = varType(0), cvX = varType(0), cvY = varType(0);
        accumulateAffineDevice(d_u, x / dx, y / dy - varType(0.5),
                               nx + 1, ny, dx, dy, kernelRadius, u, cuX, cuY);
        accumulateAffineDevice(d_v, x / dx - varType(0.5), y / dy,
                               nx, ny + 1, dx, dy, kernelRadius, v, cvX, cvY);
        d_vel_x[p] = u;
        d_vel_y[p] = v - dt * gravity;
        d_cu_x[p] = cuX;
        d_cu_y[p] = cuY;
        d_cv_x[p] = cvX;
        d_cv_y[p] = cvY;
        return;
    }

    d_vel_x[p] = uNew;
    d_vel_y[p] = vNew - dt * gravity;
}

// ─── Advection ────────────────────────────────────────────────────────────────

__global__ void advectKernel(
    int nParticles, int nx, int ny,
    varType dx, varType dy, varType dt,
    const varType* __restrict__ d_u,
    const varType* __restrict__ d_v,
    const uint16_t* __restrict__ d_labels,
    const varType* __restrict__ d_vel_x,
    const varType* __restrict__ d_vel_y,
    varType* d_pos_x, varType* d_pos_y,
    int* d_cell_idx,
    int* d_alive)
{
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= nParticles) return;

    const varType x0 = d_pos_x[p];
    const varType y0 = d_pos_y[p];
    const varType u0 = d_vel_x[p];
    const varType v0 = d_vel_y[p];

#ifdef USE_SPEED
    varType x1 = x0 + dt * u0;
    varType y1 = y0 + dt * v0;
#else
    const varType xm   = x0 + varType(0.5) * dt * u0;
    const varType ym   = y0 + varType(0.5) * dt * v0;
    const varType umid = bilinear_device<0>(xm, ym, dx, dy, d_u, nx + 1, ny);
    const varType vmid = bilinear_device<1>(xm, ym, dx, dy, d_v, nx,     ny + 1);

    varType x1 = x0 + dt * umid;
    varType y1 = y0 + dt * vmid;
#endif

    if (x1 < varType(0) || x1 >= dx * varType(nx) ||
        y1 < varType(0) || y1 >= dy * varType(ny)) {
        d_alive[p] = 0;
        return;
    }

    if (!isFluidPoint_device(x1, y1, nx, ny, dx, dy, d_labels))
        clampToLastFluidPoint_device(x0, y0, x1, y1, x1, y1,
                                     nx, ny, dx, dy, d_labels);

    d_pos_x[p]    = x1;
    d_pos_y[p]    = y1;
    d_cell_idx[p] = min((int)(x1 / dx), nx - 1)
                  + nx * min((int)(y1 / dy), ny - 1);
    d_alive[p]    = 1;
}

// ─── Cell count ───────────────────────────────────────────────────────────────
// Counts particles per cell via atomicAdd. O(N) — no sorting required.
// d_cell_start is not maintained: only d_cell_count is needed downstream.

__global__ void countCellsKernel(
    int nParticles, int nCells,
    const int* __restrict__ d_cell_idx,
    int* d_cell_count)
{
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= nParticles) return;

    const int c = d_cell_idx[p];
    if (c >= 0 && c < nCells) atomicAdd(&d_cell_count[c], 1);
}

// ─── Cell label update ────────────────────────────────────────────────────────

__global__ void cellLabelKernel(
    int nx, int ny,
    const int* __restrict__ d_cell_count,
    uint16_t* d_labels)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || j >= ny) return;

    const int lbl_idx = (nx + 2) * (j + 1) + (i + 1);
    const uint16_t cur = d_labels[lbl_idx];

    if (cur & Fields2D::SOLID) return;
    if (cur & Fields2D::BC_U)  return;
    if (cur & Fields2D::BC_V)  return;

    const uint16_t base = cur & ~uint16_t(Fields2D::SOLID | Fields2D::AIR);
    d_labels[lbl_idx] = (d_cell_count[nx * j + i] > 0)
                      ? base
                      : base | uint16_t(Fields2D::AIR);
}

// ─── GPU pressure solve: compute RHS ─────────────────────────────────────────
// b[i,j] = -coef * div(u,v) for each fluid cell. 0 for non-fluid cells.

__global__ void computeRHSKernel(
    int nx, int ny, varType dx, varType dy, varType coef,
    const uint16_t* __restrict__ d_labels,
    const varType* __restrict__ d_u,   // (nx+1)*ny
    const varType* __restrict__ d_v,   // nx*(ny+1)
    varType* d_rhs)                    // (nx+2)*(ny+2)
{
    const int i   = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const int j   = blockIdx.y * blockDim.y + threadIdx.y + 1;
    const int pnx = nx + 2;

    if (i > nx || j > ny) return;

    const int      id  = pnx * j + i;
    const uint16_t lbl = d_labels[id];

    if ((lbl & Fields2D::SOLID) || (lbl & Fields2D::AIR) || (lbl & Fields2D::BC_P)) {
        d_rhs[id] = varType(0);
        return;
    }

    const int ci = i - 1, cj = j - 1;          // 0-based physical cell

    const varType u_r = d_u[(nx + 1) * cj + (ci + 1)];
    const varType u_l = d_u[(nx + 1) * cj + ci];
    const varType v_t = d_v[nx * (cj + 1) + ci];
    const varType v_b = d_v[nx * cj + ci];

    d_rhs[id] = -coef * ((u_r - u_l) / dx + (v_t - v_b) / dy);
}

// ─── CG dot product → device ─────────────────────────────────────────────────
// Computes a·b (double precision) and atomicAdds the result into *out.
// Caller must zero *out (cudaMemset) before launching.
// Requires compute capability 6.0+ for double atomicAdd (all modern GPUs).

__global__ void dotToDevKernel(
    int n,
    const varType* __restrict__ a,
    const varType* __restrict__ b,
    double* out)
{
    __shared__ double smem[BLOCK];
    double s = 0.0;
    for (int i = blockIdx.x * BLOCK + threadIdx.x; i < n; i += gridDim.x * BLOCK)
        s += (double)a[i] * (double)b[i];
    smem[threadIdx.x] = s;
    __syncthreads();
    for (int k = BLOCK / 2; k > 0; k >>= 1) {
        if (threadIdx.x < k) smem[threadIdx.x] += smem[threadIdx.x + k];
        __syncthreads();
    }
    if (threadIdx.x == 0) atomicAdd(out, smem[0]);
}

// ─── CG pressure kernels ──────────────────────────────────────────────────────
// The pressure Laplacian stencil: (Ax)[id] = 4*x[id] - sum_fluid_neighbors(x).
// Solid/BC_U|V neighbours contribute x[id] (Neumann); AIR neighbours
// contribute 0 (Dirichlet p=0); fluid neighbours contribute x[nb].

__global__ void cgApplyLaplacianKernel(
    int pnx, int pny, int pNxy,
    const uint16_t* __restrict__ d_labels,
    const varType* __restrict__ x,
    varType* y)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= pNxy) return;

    const int i = id % pnx;
    const int j = id / pnx;

    if (i == 0 || i == pnx - 1 || j == 0 || j == pny - 1) { y[id] = varType(0); return; }

    const uint16_t lbl = d_labels[id];
    if ((lbl & Fields2D::SOLID) || (lbl & Fields2D::AIR) || (lbl & Fields2D::BC_P)) {
        y[id] = varType(0); return;
    }

    const varType xij = x[id];
    varType sumN = varType(0);  // sum of neighbour contributions

    // Left
    { const uint16_t nb = d_labels[pnx * j + (i - 1)];
      if      ((nb & Fields2D::SOLID) || (nb & Fields2D::BC_U)) sumN += xij;
      else if (!(nb & Fields2D::AIR))                           sumN += x[pnx * j + (i - 1)]; }
    // Right  (BC_U lives on current cell's right face)
    { const uint16_t nb = d_labels[pnx * j + (i + 1)];
      if      ((nb & Fields2D::SOLID) || (lbl & Fields2D::BC_U)) sumN += xij;
      else if (!(nb & Fields2D::AIR))                            sumN += x[pnx * j + (i + 1)]; }
    // Bottom
    { const uint16_t nb = d_labels[pnx * (j - 1) + i];
      if      ((nb & Fields2D::SOLID) || (nb & Fields2D::BC_V)) sumN += xij;
      else if (!(nb & Fields2D::AIR))                           sumN += x[pnx * (j - 1) + i]; }
    // Top    (BC_V lives on current cell's top face)
    { const uint16_t nb = d_labels[pnx * (j + 1) + i];
      if      ((nb & Fields2D::SOLID) || (lbl & Fields2D::BC_V)) sumN += xij;
      else if (!(nb & Fields2D::AIR))                            sumN += x[pnx * (j + 1) + i]; }

    y[id] = varType(4) * xij - sumN;
}

// r = b - Ax;  d = r.  Non-fluid cells are zeroed.
__global__ void cgInitKernel(
    int pNxy,
    const uint16_t* __restrict__ d_labels,
    const varType* __restrict__ b,
    const varType* __restrict__ Ax,
    varType* r, varType* d_dir)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= pNxy) return;

    const uint16_t lbl = d_labels[id];
    if ((lbl & Fields2D::SOLID) || (lbl & Fields2D::AIR) || (lbl & Fields2D::BC_P)) {
        r[id] = varType(0); d_dir[id] = varType(0); return;
    }
    const varType val = b[id] - Ax[id];
    r[id] = val; d_dir[id] = val;
}

// p += alpha * d;   r -= alpha * q.
// alpha = scalars[0] / scalars[1]  (sigma / dAd), computed on device.
__global__ void cgUpdatePRKernel(
    int pNxy,
    const uint16_t* __restrict__ d_labels,
    const double* __restrict__ scalars,  // [0]=sigma, [1]=dAd
    const varType* __restrict__ d_dir,
    const varType* __restrict__ q,
    varType* p, varType* r)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= pNxy) return;

    const uint16_t lbl = d_labels[id];
    if ((lbl & Fields2D::SOLID) || (lbl & Fields2D::AIR) || (lbl & Fields2D::BC_P)) return;

    const double dAd = scalars[1];
    if (!(fabs(dAd) > 1e-300)) return;
    const varType alpha = (varType)(scalars[0] / dAd);
    p[id] += alpha * d_dir[id];
    r[id] -= alpha * q[id];
}

// d = r + beta * d.
// beta = scalars[2] / scalars[0]  (sigma_new / sigma), computed on device.
__global__ void cgUpdateDirKernel(
    int pNxy,
    const uint16_t* __restrict__ d_labels,
    const double* __restrict__ scalars,  // [0]=sigma, [2]=sigma_new
    const varType* __restrict__ r,
    varType* d_dir)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= pNxy) return;

    const uint16_t lbl = d_labels[id];
    if ((lbl & Fields2D::SOLID) || (lbl & Fields2D::AIR) || (lbl & Fields2D::BC_P)) {
        d_dir[id] = varType(0); return;
    }
    const double sigma = scalars[0];
    const varType beta = (fabs(sigma) > 1e-300) ? (varType)(scalars[2] / sigma) : varType(0);
    d_dir[id] = r[id] + beta * d_dir[id];
}

// ─── Pressure gradient: U-faces ──────────────────────────────────────────────

__global__ void applyPressureGradKernel_U(
    int nx, int ny, varType dx, varType dt, varType density,
    const varType* __restrict__ d_p,
    const uint16_t* __restrict__ d_labels,
    varType* d_u)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= (nx + 1) * ny) return;

    const int fi  = id % (nx + 1);
    const int fj  = id / (nx + 1);
    const int pnx = nx + 2;

    const uint16_t lbl_L = d_labels[pnx * (fj + 1) + fi];
    const uint16_t lbl_R = d_labels[pnx * (fj + 1) + fi + 1];

    if (lbl_L & Fields2D::BC_U)                                  return;
    if ((lbl_L & Fields2D::SOLID) || (lbl_R & Fields2D::SOLID))  return;
    if ((lbl_L & Fields2D::AIR)   && (lbl_R & Fields2D::AIR))    return;

    const varType p_L = (lbl_L & Fields2D::AIR) ? varType(0) : d_p[pnx * (fj + 1) + fi];
    const varType p_R = (lbl_R & Fields2D::AIR) ? varType(0) : d_p[pnx * (fj + 1) + fi + 1];

    d_u[id] -= (dt / density) * (p_R - p_L) / dx;
}

// ─── Pressure gradient: V-faces ──────────────────────────────────────────────

__global__ void applyPressureGradKernel_V(
    int nx, int ny, varType dy, varType dt, varType density,
    const varType* __restrict__ d_p,
    const uint16_t* __restrict__ d_labels,
    varType* d_v)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= nx * (ny + 1)) return;

    const int fi  = id % nx;
    const int fj  = id / nx;
    const int pnx = nx + 2;

    const uint16_t lbl_B = d_labels[pnx * fj       + (fi + 1)];
    const uint16_t lbl_T = d_labels[pnx * (fj + 1) + (fi + 1)];

    if (lbl_B & Fields2D::BC_V)                                  return;
    if ((lbl_B & Fields2D::SOLID) || (lbl_T & Fields2D::SOLID))  return;
    if ((lbl_B & Fields2D::AIR)   && (lbl_T & Fields2D::AIR))    return;

    const varType p_B = (lbl_B & Fields2D::AIR) ? varType(0) : d_p[pnx * fj       + (fi + 1)];
    const varType p_T = (lbl_T & Fields2D::AIR) ? varType(0) : d_p[pnx * (fj + 1) + (fi + 1)];

    d_v[id] -= (dt / density) * (p_T - p_B) / dy;
}

// ─── Air cleanup ─────────────────────────────────────────────────────────────
// GPU version of the CPU air-cleanup loop.
// Zeros p at air cells; zeros u/v at faces adjacent to non-fluid cells.
// Multiple threads may write 0 to the same face (idempotent, benign race).

__global__ void airCleanupKernel(
    int nx, int ny,
    const uint16_t* __restrict__ d_labels,
    varType* d_u, varType* d_v, varType* d_p)
{
    const int i   = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const int j   = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if (i > nx || j > ny) return;

    const int      pnx = nx + 2;
    const int      id  = pnx * j + i;
    const uint16_t lbl = d_labels[id];

    if (!(lbl & Fields2D::AIR)) return;

    d_p[id] = varType(0);

    const int ci = i - 1, cj = j - 1;

    const uint16_t lL = d_labels[pnx * j       + (i - 1)];
    const uint16_t lR = d_labels[pnx * j       + (i + 1)];
    const uint16_t lB = d_labels[pnx * (j - 1) + i];
    const uint16_t lT = d_labels[pnx * (j + 1) + i];

    if (!IS_FLUID(lL)) d_u[(nx + 1) * cj + ci]       = varType(0);
    if (!IS_FLUID(lR)) d_u[(nx + 1) * cj + (ci + 1)] = varType(0);
    if (!IS_FLUID(lB)) d_v[nx * cj + ci]              = varType(0);
    if (!IS_FLUID(lT)) d_v[nx * (cj + 1) + ci]        = varType(0);
}

// ─── GPU refill ──────────────────────────────────────────────────────────────

__device__ __forceinline__ uint32_t hashU32(uint32_t x) {
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

__device__ __forceinline__ varType rand01_device(int cell, int particle, int seed, int salt) {
    const uint32_t h = hashU32(static_cast<uint32_t>(cell) * 0x9e3779b9u
                             ^ static_cast<uint32_t>(particle) * 0x85ebca6bu
                             ^ static_cast<uint32_t>(seed) * 0xc2b2ae35u
                             ^ static_cast<uint32_t>(salt));
    return static_cast<varType>((h >> 8) * (1.0 / 16777216.0));
}

__global__ void refillKernel(
    int nx, int ny, int targetPPC, int maxParticles, int seed,
    varType dx, varType dy,
    const uint16_t* __restrict__ d_labels,
    const int* __restrict__ d_cell_count,
    const varType* __restrict__ d_u,
    const varType* __restrict__ d_v,
    varType* d_pos_x,
    varType* d_pos_y,
    varType* d_vel_x,
    varType* d_vel_y,
    varType* d_cu_x,
    varType* d_cu_y,
    varType* d_cv_x,
    varType* d_cv_y,
    int* d_cell_idx,
    int* d_alive,
    int* d_particle_count,
    int* d_refill_overflow)
{
    const int c = blockIdx.x * blockDim.x + threadIdx.x;
    const int nCells = nx * ny;
    if (c >= nCells) return;

    const int ci = c % nx;
    const int cj = c / nx;
    const uint16_t lbl = d_labels[(nx + 2) * (cj + 1) + (ci + 1)];
    if ((lbl & Fields2D::SOLID) || (lbl & Fields2D::AIR)) return;

    const int missing = targetPPC - d_cell_count[c];
    if (missing <= 0) return;

    const int base = atomicAdd(d_particle_count, missing);
    const int available = maxParticles - base;
    if (available <= 0) {
        atomicAdd(d_refill_overflow, missing);
        return;
    }

    const int toWrite = min(missing, available);
    if (toWrite < missing)
        atomicAdd(d_refill_overflow, missing - toWrite);

    for (int m = 0; m < toWrite; ++m) {
        const int slot = base + m;
        const varType rx = rand01_device(c, m, seed, 0x1234u);
        const varType ry = rand01_device(c, m, seed, 0x9abcu);
        const varType x = (static_cast<varType>(ci) + rx) * dx;
        const varType y = (static_cast<varType>(cj) + ry) * dy;

        d_pos_x[slot] = x;
        d_pos_y[slot] = y;
        d_vel_x[slot] = bilinear_device<0>(x, y, dx, dy, d_u, nx + 1, ny);
        d_vel_y[slot] = bilinear_device<1>(x, y, dx, dy, d_v, nx, ny + 1);
        if (d_cu_x) {
            d_cu_x[slot] = varType(0);
            d_cu_y[slot] = varType(0);
            d_cv_x[slot] = varType(0);
            d_cv_y[slot] = varType(0);
        }
        d_cell_idx[slot] = c;
        d_alive[slot] = 1;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// GPIC host implementation
// ═══════════════════════════════════════════════════════════════════════════════

// ─── Constructor ──────────────────────────────────────────────────────────────

GPIC::GPIC(Parameters& params)
    : Solver(params),
      dev_(std::make_unique<DeviceState>())
{
    // Build Cloud2D as VanillaPIC so needsAffine == false.
    const auto savedMethod = params.solver.method;
    params.solver.method   = SolverConfig::Method::VanillaPIC;
    cloud_ = std::make_unique<Cloud2D>(params);
    params.solver.method   = savedMethod;
    params.applyToFields(*fields);
    cloud_->InitParticleGrid(*fields, params.ppcx, params.ppcy);

    if (params.write_particles)
        particlesWriter_ = std::make_unique<OutputWriter>(params.folder, "particles");

    auto& d    = *dev_;
    d.nx       = nx;  d.ny = ny;
    d.dx       = dx;  d.dy = dy;  d.dt = dt;
    d.density  = static_cast<varType>(params.density);
    if (params.solver.method == SolverConfig::Method::GFLIP) {
        d.transferMode = TRANSFER_FLIP;
    } else if (params.solver.method == SolverConfig::Method::GAPIC) {
        d.transferMode = TRANSFER_APIC;
        d.needsAffine = true;
    }

    const int total    = cloud_->totalSize();
    d.maxParticles     = std::max(total * 4, 64);
    d.nParticles       = total;

    const std::size_t nP   = static_cast<std::size_t>(d.maxParticles);
    const std::size_t u_sz = static_cast<std::size_t>(nx + 1) * ny;
    const std::size_t v_sz = static_cast<std::size_t>(nx) * (ny + 1);
    const std::size_t p_sz = static_cast<std::size_t>(nx + 2) * (ny + 2);
    const std::size_t c_sz = static_cast<std::size_t>(nx) * ny;

    // ── Particle SoA ──
    CUDA_CHECK(cudaMalloc(&d.d_pos_x,    nP * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_pos_y,    nP * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_vel_x,    nP * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_vel_y,    nP * sizeof(varType)));
    if (d.needsAffine) {
        CUDA_CHECK(cudaMalloc(&d.d_cu_x, nP * sizeof(varType)));
        CUDA_CHECK(cudaMalloc(&d.d_cu_y, nP * sizeof(varType)));
        CUDA_CHECK(cudaMalloc(&d.d_cv_x, nP * sizeof(varType)));
        CUDA_CHECK(cudaMalloc(&d.d_cv_y, nP * sizeof(varType)));
        CUDA_CHECK(cudaMemset(d.d_cu_x, 0, nP * sizeof(varType)));
        CUDA_CHECK(cudaMemset(d.d_cu_y, 0, nP * sizeof(varType)));
        CUDA_CHECK(cudaMemset(d.d_cv_x, 0, nP * sizeof(varType)));
        CUDA_CHECK(cudaMemset(d.d_cv_y, 0, nP * sizeof(varType)));
    }
    CUDA_CHECK(cudaMalloc(&d.d_cell_idx, nP * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d.d_alive,    nP * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d.d_tmp,      nP * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_perm,     nP * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d.d_particle_count,  sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d.d_refill_overflow, sizeof(int)));

    // ── Grid velocity ──
    CUDA_CHECK(cudaMalloc(&d.d_u,     u_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_v,     v_sz * sizeof(varType)));
    if (d.transferMode == TRANSFER_FLIP) {
        CUDA_CHECK(cudaMalloc(&d.d_u_old, u_sz * sizeof(varType)));
        CUDA_CHECK(cudaMalloc(&d.d_v_old, v_sz * sizeof(varType)));
    }
    CUDA_CHECK(cudaMalloc(&d.d_u_num, u_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_u_den, u_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_v_num, v_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_v_den, v_sz * sizeof(varType)));

    // ── Pressure / CG workspace ──
    CUDA_CHECK(cudaMalloc(&d.d_p,          p_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_rhs,        p_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_cg_r,       p_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_cg_d,       p_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_cg_q,       p_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_cg_scalars, 3 * sizeof(double)));
    CUDA_CHECK(cudaMemset(d.d_p,           0, p_sz * sizeof(varType)));

    // ── Labels + cell occupancy ──
    CUDA_CHECK(cudaMalloc(&d.d_labels,     p_sz * sizeof(uint16_t)));
    CUDA_CHECK(cudaMalloc(&d.d_cell_count, c_sz * sizeof(int)));
}

// ─── Destructor ───────────────────────────────────────────────────────────────

GPIC::~GPIC() {
    if (!dev_) return;
    auto& d = *dev_;
#define FREE(p) if (p) { cudaFree(p); (p) = nullptr; }
    FREE(d.d_pos_x)     FREE(d.d_pos_y)
    FREE(d.d_vel_x)     FREE(d.d_vel_y)
    FREE(d.d_cu_x)      FREE(d.d_cu_y)
    FREE(d.d_cv_x)      FREE(d.d_cv_y)
    FREE(d.d_cell_idx)  FREE(d.d_alive)
    FREE(d.d_tmp)       FREE(d.d_perm)
    FREE(d.d_particle_count) FREE(d.d_refill_overflow)
    FREE(d.d_u)         FREE(d.d_v)
    FREE(d.d_u_old)     FREE(d.d_v_old)
    FREE(d.d_u_num)     FREE(d.d_u_den)
    FREE(d.d_v_num)     FREE(d.d_v_den)
    FREE(d.d_p)         FREE(d.d_rhs)
    FREE(d.d_cg_r)      FREE(d.d_cg_d)      FREE(d.d_cg_q)
    FREE(d.d_cg_scalars)
    FREE(d.d_labels)    FREE(d.d_cell_count)
#undef FREE
}

// ─── Upload helpers ───────────────────────────────────────────────────────────

void GPIC::UploadParticles() {
    auto& d = *dev_;

    std::vector<varType> h_px, h_py, h_vx, h_vy;
    std::vector<int>     h_ci;

    for (int cj = 0; cj < ny; ++cj) {
        for (int ci = 0; ci < nx; ++ci) {
            const Particles& cell     = (*cloud_)(ci, cj);
            const int        nc       = cell.size();
            const int        cell_idx = nx * cj + ci;
            for (int p = 0; p < nc; ++p) {
                h_px.push_back(cell.GetX(p));
                h_py.push_back(cell.GetY(p));
                h_vx.push_back(cell.GetU(p));
                h_vy.push_back(cell.GetV(p));
                h_ci.push_back(cell_idx);
            }
        }
    }

    const int N = static_cast<int>(h_px.size());
    d.nParticles = N;
    if (N == 0) return;

    if (N > d.maxParticles) {
        std::fprintf(stderr, "[GPIC] particle count %d > maxParticles %d\n", N, d.maxParticles);
        std::exit(EXIT_FAILURE);
    }

    CUDA_CHECK(cudaMemcpy(d.d_pos_x,    h_px.data(), N * sizeof(varType), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.d_pos_y,    h_py.data(), N * sizeof(varType), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.d_vel_x,    h_vx.data(), N * sizeof(varType), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.d_vel_y,    h_vy.data(), N * sizeof(varType), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.d_cell_idx, h_ci.data(), N * sizeof(int),     cudaMemcpyHostToDevice));
    if (d.needsAffine) {
        CUDA_CHECK(cudaMemset(d.d_cu_x, 0, N * sizeof(varType)));
        CUDA_CHECK(cudaMemset(d.d_cu_y, 0, N * sizeof(varType)));
        CUDA_CHECK(cudaMemset(d.d_cv_x, 0, N * sizeof(varType)));
        CUDA_CHECK(cudaMemset(d.d_cv_y, 0, N * sizeof(varType)));
    }
}

void GPIC::UploadGrid() {
    auto& d = *dev_;
    const std::size_t u_sz = static_cast<std::size_t>(nx + 1) * ny;
    const std::size_t v_sz = static_cast<std::size_t>(nx) * (ny + 1);
    const std::size_t p_sz = static_cast<std::size_t>(nx + 2) * (ny + 2);

    CUDA_CHECK(cudaMemcpy(d.d_u,      fields->u.A,    u_sz * sizeof(varType),  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.d_v,      fields->v.A,    v_sz * sizeof(varType),  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.d_labels, fields->Labels.get(), p_sz * sizeof(uint16_t), cudaMemcpyHostToDevice));
}

// ─── Download helpers ─────────────────────────────────────────────────────────

void GPIC::DownloadParticles() const {
    auto& d = *dev_;
    const int N = d.nParticles;

    for (auto& cell : cloud_->cells) {
        cell.pos_x.clear(); cell.pos_y.clear();
        cell.vel_x.clear(); cell.vel_y.clear();
        if (cell.needsAffine) {
            cell.cu_x.clear(); cell.cu_y.clear();
            cell.cv_x.clear(); cell.cv_y.clear();
        }
    }

    if (N == 0) return;

    std::vector<varType> h_px(N), h_py(N), h_vx(N), h_vy(N);
    std::vector<int>     h_ci(N);

    CUDA_CHECK(cudaMemcpy(h_px.data(), d.d_pos_x,    N * sizeof(varType), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_py.data(), d.d_pos_y,    N * sizeof(varType), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vx.data(), d.d_vel_x,    N * sizeof(varType), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vy.data(), d.d_vel_y,    N * sizeof(varType), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_ci.data(), d.d_cell_idx, N * sizeof(int),     cudaMemcpyDeviceToHost));

    for (int p = 0; p < N; ++p) {
        const int c = h_ci[p];
        if (c < 0 || c >= nx * ny) continue;
        (*cloud_)(c % nx, c / nx).Add(h_px[p], h_py[p], h_vx[p], h_vy[p], 0u);
    }
}

void GPIC::DownloadAll() const {
    auto& d = *dev_;
    const std::size_t u_sz = static_cast<std::size_t>(nx + 1) * ny;
    const std::size_t v_sz = static_cast<std::size_t>(nx) * (ny + 1);
    const std::size_t p_sz = static_cast<std::size_t>(nx + 2) * (ny + 2);

    CUDA_CHECK(cudaMemcpy(fields->u.A,    d.d_u,      u_sz * sizeof(varType),  cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(fields->v.A,    d.d_v,      v_sz * sizeof(varType),  cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(fields->p.A,    d.d_p,      p_sz * sizeof(varType),  cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(fields->Labels.get(), d.d_labels, p_sz * sizeof(uint16_t), cudaMemcpyDeviceToHost));
}

void GPIC::DownloadLabels() const {
    auto& d = *dev_;
    const std::size_t p_sz = static_cast<std::size_t>(nx + 2) * (ny + 2);
    CUDA_CHECK(cudaMemcpy(fields->Labels.get(), d.d_labels,
                          p_sz * sizeof(uint16_t), cudaMemcpyDeviceToHost));
}

// ─── GPU step stages ──────────────────────────────────────────────────────────

void GPIC::GPUProjectParticlesOnGrid() {
    auto& d = *dev_;
    const int N = d.nParticles;
    const std::size_t u_sz = static_cast<std::size_t>(nx + 1) * ny;
    const std::size_t v_sz = static_cast<std::size_t>(nx) * (ny + 1);

    CUDA_CHECK(cudaMemset(d.d_u_num, 0, u_sz * sizeof(varType)));
    CUDA_CHECK(cudaMemset(d.d_u_den, 0, u_sz * sizeof(varType)));
    CUDA_CHECK(cudaMemset(d.d_v_num, 0, v_sz * sizeof(varType)));
    CUDA_CHECK(cudaMemset(d.d_v_den, 0, v_sz * sizeof(varType)));

    if (N > 0) {
        const int g = (N + BLOCK - 1) / BLOCK;
        p2gKernel<<<g, BLOCK>>>(
            N, nx, ny, d.transferMode, dx, dy,
            d.d_pos_x, d.d_pos_y, d.d_vel_x, d.d_vel_y,
            d.d_cu_x, d.d_cu_y, d.d_cv_x, d.d_cv_y,
            d.d_u_num, d.d_u_den, d.d_v_num, d.d_v_den);
        CUDA_CHECK(cudaGetLastError());
    }

    {
        const int gu = ((nx + 1) * ny + BLOCK - 1) / BLOCK;
        p2gNormaliseKernel_U<<<gu, BLOCK>>>(
            nx, ny, d.d_u_num, d.d_u_den, d.d_labels, d.d_u);
        CUDA_CHECK(cudaGetLastError());

        const int gv = (nx * (ny + 1) + BLOCK - 1) / BLOCK;
        p2gNormaliseKernel_V<<<gv, BLOCK>>>(
            nx, ny, d.d_v_num, d.d_v_den, d.d_labels, d.d_v);
        CUDA_CHECK(cudaGetLastError());
    }
}

// ─── GPU pressure solve (CG) ──────────────────────────────────────────────────
// 1. Compute RHS b = -coef * div(u,v) on GPU.
// 2. Solve A*p = b with unpreconditioned CG (warm start from previous step).
// 3. Apply pressure gradient to d_u, d_v on GPU.
//
// CG converges in O(sqrt(N)) iterations for the pressure Laplacian vs O(N)
// for RBGS, giving a large speedup on grids that are not too large.
// All three CG vectors (r, d, q) live in persistent device memory and are
// re-initialised each step from the warm-start pressure.

void GPIC::GPUMakeIncompressible() {
    auto& d = *dev_;
    const int pnx  = nx + 2;
    const int pny  = ny + 2;
    const int pNxy = pnx * pny;
    const varType coef = d.density * dx * dx / d.dt;

    // ── 1. Compute RHS b ──
    CUDA_CHECK(cudaMemset(d.d_rhs, 0, pNxy * sizeof(varType)));
    {
        const dim3 block(BX, BY);
        const dim3 grid((nx + BX - 1) / BX, (ny + BY - 1) / BY);
        computeRHSKernel<<<grid, block>>>(
            nx, ny, dx, dy, coef, d.d_labels, d.d_u, d.d_v, d.d_rhs);
        CUDA_CHECK(cudaGetLastError());
    }

    // ── 2. ||b||² — one D2H sync, early exit when already divergence-free ──
    {
        auto rhs_ptr = thrust::device_pointer_cast(d.d_rhs);
        const varType bNormSq = thrust::transform_reduce(
            rhs_ptr, rhs_ptr + pNxy, SquareFunctor{}, varType(0),
            cuda::std::plus<varType>{});

        if (bNormSq >= REAL_EPSILON * REAL_EPSILON) {
            const double tolSq = (double)params.solver.tolerance
                               * (double)params.solver.tolerance
                               * (double)bNormSq;

            // ── 3. CG solve (warm start on d_p) ──
            // All dot products write to d_cg_scalars on device.
            // D2H sync only every CG_CHECK iterations for convergence.
            const int cgBlocks = (pNxy + BLOCK - 1) / BLOCK;
            const int maxIters = params.solver.maxIters;

            // q = A * p  (initial)
            cgApplyLaplacianKernel<<<cgBlocks, BLOCK>>>(pnx, pny, pNxy, d.d_labels, d.d_p, d.d_cg_q);
            CUDA_CHECK(cudaGetLastError());

            // r = b - A*p;   d_dir = r
            cgInitKernel<<<cgBlocks, BLOCK>>>(pNxy, d.d_labels, d.d_rhs, d.d_cg_q, d.d_cg_r, d.d_cg_d);
            CUDA_CHECK(cudaGetLastError());

            // scalars[0] = sigma = r · r  (device)
            CUDA_CHECK(cudaMemset(d.d_cg_scalars, 0, sizeof(double)));
            dotToDevKernel<<<DOT_GRIDS, BLOCK>>>(pNxy, d.d_cg_r, d.d_cg_r, d.d_cg_scalars);

            for (int it = 0; it < maxIters; ++it) {
                // Convergence check: D2H every CG_CHECK iters (scalars[0] = current sigma)
                if (it % CG_CHECK == 0) {
                    double sigma_h;
                    CUDA_CHECK(cudaMemcpy(&sigma_h, d.d_cg_scalars, sizeof(double),
                                         cudaMemcpyDeviceToHost));
                    if (!std::isfinite(sigma_h) || sigma_h <= tolSq) break;
                }

                // q = A * d
                cgApplyLaplacianKernel<<<cgBlocks, BLOCK>>>(pnx, pny, pNxy,
                    d.d_labels, d.d_cg_d, d.d_cg_q);
                CUDA_CHECK(cudaGetLastError());

                // scalars[1] = dAd = d · q  (device)
                CUDA_CHECK(cudaMemset(d.d_cg_scalars + 1, 0, sizeof(double)));
                dotToDevKernel<<<DOT_GRIDS, BLOCK>>>(pNxy, d.d_cg_d, d.d_cg_q, d.d_cg_scalars + 1);

                // p += (sigma/dAd)*d;   r -= (sigma/dAd)*q  — alpha computed on device
                cgUpdatePRKernel<<<cgBlocks, BLOCK>>>(pNxy, d.d_labels, d.d_cg_scalars,
                    d.d_cg_d, d.d_cg_q, d.d_p, d.d_cg_r);
                CUDA_CHECK(cudaGetLastError());

                // scalars[2] = sigma_new = r · r  (device)
                CUDA_CHECK(cudaMemset(d.d_cg_scalars + 2, 0, sizeof(double)));
                dotToDevKernel<<<DOT_GRIDS, BLOCK>>>(pNxy, d.d_cg_r, d.d_cg_r, d.d_cg_scalars + 2);

                // d = r + (sigma_new/sigma)*d  — beta computed on device
                cgUpdateDirKernel<<<cgBlocks, BLOCK>>>(pNxy, d.d_labels, d.d_cg_scalars,
                    d.d_cg_r, d.d_cg_d);
                CUDA_CHECK(cudaGetLastError());

                // Shift: scalars[0] = scalars[2]  (D2D copy, no PCIe)
                CUDA_CHECK(cudaMemcpy(d.d_cg_scalars, d.d_cg_scalars + 2, sizeof(double),
                                      cudaMemcpyDeviceToDevice));
            }
        }
    }

    // ── 4. Apply pressure gradient ──
    {
        const int gu = ((nx + 1) * ny + BLOCK - 1) / BLOCK;
        applyPressureGradKernel_U<<<gu, BLOCK>>>(
            nx, ny, dx, d.dt, d.density, d.d_p, d.d_labels, d.d_u);
        CUDA_CHECK(cudaGetLastError());

        const int gv = (nx * (ny + 1) + BLOCK - 1) / BLOCK;
        applyPressureGradKernel_V<<<gv, BLOCK>>>(
            nx, ny, dy, d.dt, d.density, d.d_p, d.d_labels, d.d_v);
        CUDA_CHECK(cudaGetLastError());
    }
}

void GPIC::GPUProjectGridOnParticles() {
    auto& d = *dev_;
    const int N = d.nParticles;
    if (N == 0) return;

    const int g = (N + BLOCK - 1) / BLOCK;
    g2pKernel<<<g, BLOCK>>>(
        N, nx, ny, d.transferMode, params.kernelOrder,
        dx, dy, d.dt, params.gravity, params.coefPic,
        d.d_u, d.d_v, d.d_u_old, d.d_v_old,
        d.d_pos_x, d.d_pos_y,
        d.d_vel_x, d.d_vel_y,
        d.d_cu_x, d.d_cu_y, d.d_cv_x, d.d_cv_y);
    CUDA_CHECK(cudaGetLastError());
}

void GPIC::GPUAdvect() {
    auto& d = *dev_;
    const int N = d.nParticles;
    if (N == 0) return;

    {
        const int g = (N + BLOCK - 1) / BLOCK;
        advectKernel<<<g, BLOCK>>>(
            N, nx, ny, dx, dy, d.dt,
            d.d_u, d.d_v, d.d_labels,
            d.d_vel_x, d.d_vel_y,
            d.d_pos_x, d.d_pos_y,
            d.d_cell_idx, d.d_alive);
        CUDA_CHECK(cudaGetLastError());
    }

    auto alive_ptr = thrust::device_pointer_cast(d.d_alive);
    auto tmp_ptr   = thrust::device_pointer_cast(d.d_tmp);

    const int newN = static_cast<int>(
        thrust::count(thrust::device, alive_ptr, alive_ptr + N, 1));

    if (newN < N) {
        auto copy_compact_vt = [&](varType* src) {
            thrust::copy_if(thrust::device,
                thrust::device_pointer_cast(src),
                thrust::device_pointer_cast(src + N),
                alive_ptr, tmp_ptr, NonZero{});
            CUDA_CHECK(cudaMemcpy(src, d.d_tmp, newN * sizeof(varType),
                                  cudaMemcpyDeviceToDevice));
        };
        copy_compact_vt(d.d_pos_x);
        copy_compact_vt(d.d_pos_y);
        copy_compact_vt(d.d_vel_x);
        copy_compact_vt(d.d_vel_y);
        if (d.needsAffine) {
            copy_compact_vt(d.d_cu_x);
            copy_compact_vt(d.d_cu_y);
            copy_compact_vt(d.d_cv_x);
            copy_compact_vt(d.d_cv_y);
        }

        thrust::copy_if(thrust::device,
            thrust::device_pointer_cast(d.d_cell_idx),
            thrust::device_pointer_cast(d.d_cell_idx + N),
            alive_ptr,
            thrust::device_pointer_cast(d.d_perm), NonZero{});
        CUDA_CHECK(cudaMemcpy(d.d_cell_idx, d.d_perm, newN * sizeof(int),
                              cudaMemcpyDeviceToDevice));
    }

    d.nParticles = newN;
}

void GPIC::GPUCountParticles() {
    auto& d = *dev_;
    const int N = d.nParticles;

    CUDA_CHECK(cudaMemset(d.d_cell_count, 0, static_cast<std::size_t>(nx) * ny * sizeof(int)));

    if (N > 0) {
        const int g = (N + BLOCK - 1) / BLOCK;
        countCellsKernel<<<g, BLOCK>>>(N, nx * ny, d.d_cell_idx, d.d_cell_count);
        CUDA_CHECK(cudaGetLastError());
    }

}

// Fully on GPU: labels updated, then air cells cleared.
// No CPU download — labels stay on device until WriteOutput.
void GPIC::GPUUpdateCellState() {
    auto& d = *dev_;

    const dim3 block2d(BX, BY);
    const dim3 grid2d((nx + BX - 1) / BX, (ny + BY - 1) / BY);

    cellLabelKernel<<<grid2d, block2d>>>(nx, ny, d.d_cell_count, d.d_labels);
    CUDA_CHECK(cudaGetLastError());

    airCleanupKernel<<<grid2d, block2d>>>(nx, ny, d.d_labels, d.d_u, d.d_v, d.d_p);
    CUDA_CHECK(cudaGetLastError());
}

void GPIC::GPURefillParticles() {
    auto& d = *dev_;
    const int targetPPC = params.ppcx * params.ppcy;
    if (targetPPC <= 0 || d.nParticles >= d.maxParticles)
        return;

    CUDA_CHECK(cudaMemcpy(d.d_particle_count, &d.nParticles, sizeof(int),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d.d_refill_overflow, 0, sizeof(int)));

    const int nCells = nx * ny;
    const int g = (nCells + BLOCK - 1) / BLOCK;
    refillKernel<<<g, BLOCK>>>(
        nx, ny, targetPPC, d.maxParticles, d.refillSeed++,
        dx, dy,
        d.d_labels, d.d_cell_count,
        d.d_u, d.d_v,
        d.d_pos_x, d.d_pos_y,
        d.d_vel_x, d.d_vel_y,
        d.d_cu_x, d.d_cu_y, d.d_cv_x, d.d_cv_y,
        d.d_cell_idx, d.d_alive,
        d.d_particle_count, d.d_refill_overflow);
    CUDA_CHECK(cudaGetLastError());

    int newCount = 0;
    int overflow = 0;
    CUDA_CHECK(cudaMemcpy(&newCount, d.d_particle_count, sizeof(int),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(&overflow, d.d_refill_overflow, sizeof(int),
                          cudaMemcpyDeviceToHost));

    if (overflow > 0 || newCount > d.maxParticles) {
        std::fprintf(stderr,
                     "[GPIC] refill exceeded maxParticles (%d > %d, skipped %d). "
                     "Increase the particle capacity growth factor.\n",
                     newCount, d.maxParticles, overflow);
        std::exit(EXIT_FAILURE);
    }

    d.nParticles = newCount;
}

// ─── CFL (GPU) ────────────────────────────────────────────────────────────────

int GPIC::computeAdvectionSubstepsGPU() const {
    auto& d = *dev_;
    const int N = d.nParticles;
    if (N == 0 || params.max_cfl <= REAL_EPSILON || dt <= REAL_EPSILON) return 1;

    auto vx_ptr = thrust::device_pointer_cast(d.d_vel_x);
    auto vy_ptr = thrust::device_pointer_cast(d.d_vel_y);

    const varType max_vx = thrust::transform_reduce(
        vx_ptr, vx_ptr + N, AbsFunctor{}, varType(0), cuda::maximum<varType>{});
    const varType max_vy = thrust::transform_reduce(
        vy_ptr, vy_ptr + N, AbsFunctor{}, varType(0), cuda::maximum<varType>{});

    const varType max_vy_eff = max_vy + static_cast<varType>(std::abs(params.gravity)) * dt;
    const varType maxCourant = std::max(max_vx / dx, max_vy_eff / dy);

    if (maxCourant <= REAL_EPSILON) return 1;
    return std::max(1, static_cast<int>(std::ceil(dt * maxCourant / params.max_cfl)));
}

// ─── Step ─────────────────────────────────────────────────────────────────────
// Fully on GPU. Particles and grid stay on device between substeps.
// CPU is only involved when params.refill is true (optional).

void GPIC::Step() {
    const varType frameDt = dt;
    const int substeps    = computeAdvectionSubstepsGPU();
    const varType subDt   = frameDt / static_cast<varType>(substeps);
    setTimeStep(subDt);
    dev_->dt = subDt;

    for (int s = 0; s < substeps; ++s) {
        GPUProjectParticlesOnGrid();
        if (dev_->transferMode == TRANSFER_FLIP) {
            const std::size_t u_sz = static_cast<std::size_t>(nx + 1) * ny;
            const std::size_t v_sz = static_cast<std::size_t>(nx) * (ny + 1);
            CUDA_CHECK(cudaMemcpy(dev_->d_u_old, dev_->d_u,
                                  u_sz * sizeof(varType),
                                  cudaMemcpyDeviceToDevice));
            CUDA_CHECK(cudaMemcpy(dev_->d_v_old, dev_->d_v,
                                  v_sz * sizeof(varType),
                                  cudaMemcpyDeviceToDevice));
        }
        GPUMakeIncompressible();
        GPUProjectGridOnParticles();
        GPUAdvect();
        GPUCountParticles();
        GPUUpdateCellState();

        if (params.refill)
            GPURefillParticles();
    }

    setTimeStep(frameDt);
    dev_->dt = frameDt;
}

// ─── Run ──────────────────────────────────────────────────────────────────────

void GPIC::Run() {
    // Upload initial state once — particles and grid stay on GPU throughout
    UploadParticles();
    UploadGrid();

    const int savedTransferMode = dev_->transferMode;
    if (savedTransferMode == TRANSFER_FLIP)
        dev_->transferMode = TRANSFER_PIC;
    GPUProjectGridOnParticles();   // init particle velocities from initial grid
    dev_->transferMode = savedTransferMode;

    WriteOutput(0);                // downloads from GPU internally
    RunLoop(std::max(1, params.nt / 100));
}

// ─── WriteOutput ──────────────────────────────────────────────────────────────
// Downloads grid from GPU to CPU before calling base-class writers.

void GPIC::WriteOutput(int step) const {
    if (step % params.sampling_rate != 0) return;

    const bool needsGrid =
        params.write_u || params.write_v || params.write_p ||
        params.write_div || params.write_norm_velocity ||
        params.write_vorticity;

    if (needsGrid) {
        DownloadAll();                      // d_u, d_v, d_p, d_labels → fields
        if (params.write_div)
            fields->Div();
        if (params.write_norm_velocity)
            fields->VelocityNormCenterGrid();
        if (params.write_vorticity)
            fields->VorticityCenterGrid();
    } else {
        DownloadLabels();                   // Solver::WriteOutput always writes labels
    }

    Solver::WriteOutput(step);

    if (params.write_particles && particlesWriter_) {
        DownloadParticles();
        const bool ok = particlesWriter_->writeCloud(*cloud_, "particles");
        if (!ok)
            std::cerr << "[GPIC] Warning: failed to write particles at step "
                      << step << '\n';
    }
}

#endif // USE_CUDA
