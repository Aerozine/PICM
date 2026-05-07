#ifdef USE_CUDA

#include "GPIC.hpp"
#include "GPICKernelUtils.cuh"

#include <cuda_runtime.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/gather.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
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
    int*     d_cell_idx = nullptr;
    int*     d_alive    = nullptr;

    // Grid velocity fields (persistent on device between substeps)
    varType*  d_u      = nullptr;  // (nx+1) * ny
    varType*  d_v      = nullptr;  // nx * (ny+1)

    // Pressure + RBGS workspace (persistent = warm start between steps)
    varType*  d_p      = nullptr;  // (nx+2) * (ny+2)
    varType*  d_rhs    = nullptr;  // (nx+2) * (ny+2)  RHS for pressure solve
    varType*  d_resSq  = nullptr;  // (nx+2) * (ny+2)  per-cell residual^2

    // P2G accumulators (cleared each substep)
    varType*  d_u_num  = nullptr;  // (nx+1) * ny
    varType*  d_u_den  = nullptr;
    varType*  d_v_num  = nullptr;  // nx * (ny+1)
    varType*  d_v_den  = nullptr;

    // Labels (persistent; updated by cellLabelKernel + airCleanupKernel)
    uint16_t* d_labels = nullptr;  // (nx+2) * (ny+2)

    // Cell occupancy (rebuilt after sort)
    int* d_cell_start  = nullptr;  // [nx*ny]
    int* d_cell_count  = nullptr;  // [nx*ny]

    // Reusable scratch (avoid per-call alloc)
    varType* d_tmp     = nullptr;  // [maxParticles] varType scratch
    int*     d_perm    = nullptr;  // [maxParticles] sort permutation

    int maxParticles = 0;
    int nParticles   = 0;
    int nx = 0, ny = 0;
    varType dx{}, dy{}, dt{};
    varType density{};
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

static constexpr int BLOCK      = 256;
static constexpr int BX         = 32;
static constexpr int BY         = 8;
static constexpr int CHECK_EVERY = 32;

// ─── P2G scatter ──────────────────────────────────────────────────────────────

__global__ void p2gKernel(
    int nParticles, int nx, int ny,
    varType dx, varType dy,
    const varType* __restrict__ d_pos_x,
    const varType* __restrict__ d_pos_y,
    const varType* __restrict__ d_vel_x,
    const varType* __restrict__ d_vel_y,
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
                atomicAdd(&d_u_num[id], w * vx);
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
                atomicAdd(&d_v_num[id], w * vy);
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

__global__ void g2pKernel(
    int nParticles, int nx, int ny,
    varType dx, varType dy, varType dt, varType gravity,
    const varType* __restrict__ d_u,
    const varType* __restrict__ d_v,
    const varType* __restrict__ d_pos_x,
    const varType* __restrict__ d_pos_y,
    varType* d_vel_x,
    varType* d_vel_y)
{
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= nParticles) return;

    const varType x = d_pos_x[p];
    const varType y = d_pos_y[p];

    d_vel_x[p] = bilinear_device<0>(x, y, dx, dy, d_u, nx + 1, ny);
    d_vel_y[p] = bilinear_device<1>(x, y, dx, dy, d_v, nx,     ny + 1)
               - dt * gravity;
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

// ─── Cell ranges ─────────────────────────────────────────────────────────────

__global__ void buildCellRangesKernel(
    int nParticles, int nCells,
    const int* __restrict__ d_cell_idx,
    int* d_cell_start,
    int* d_cell_count)
{
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= nParticles) return;

    const int c = d_cell_idx[p];
    if (c < 0 || c >= nCells) return;

    atomicAdd(&d_cell_count[c], 1);
    if (p == 0 || d_cell_idx[p - 1] != c)
        d_cell_start[c] = p;
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

// ─── GPU pressure solve: RBGS sweep ──────────────────────────────────────────
// Processes one checkerboard colour per call.
// Writes squared true-residual (before update) to d_resSq for convergence.

__global__ void gpicRBGSKernel(
    int colour, int pnx, int pny, varType omega,
    const uint16_t* __restrict__ d_labels,
    const varType* __restrict__ d_rhs,
    varType* d_p,
    varType* d_resSq)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if (i >= pnx - 1 || j >= pny - 1) return;
    if (((i + j) & 1) != colour) return;

    const int      id  = pnx * j + i;
    const uint16_t lbl = d_labels[id];

    if ((lbl & Fields2D::SOLID) || (lbl & Fields2D::AIR) || (lbl & Fields2D::BC_P)) {
        d_resSq[id] = varType(0);
        return;
    }

    const varType pij = d_p[id];
    varType sumP = varType(0);

    // Left (i-1, j)
    {
        const uint16_t nb = d_labels[pnx * j + (i - 1)];
        if ((nb & Fields2D::SOLID) || (nb & Fields2D::BC_U)) sumP += pij;
        else if (!(nb & Fields2D::AIR)) sumP += d_p[pnx * j + (i - 1)];
    }
    // Right (i+1, j): BC_U checked on current cell (face belongs to cur)
    {
        const uint16_t nb = d_labels[pnx * j + (i + 1)];
        if ((nb & Fields2D::SOLID) || (lbl & Fields2D::BC_U)) sumP += pij;
        else if (!(nb & Fields2D::AIR)) sumP += d_p[pnx * j + (i + 1)];
    }
    // Bottom (i, j-1)
    {
        const uint16_t nb = d_labels[pnx * (j - 1) + i];
        if ((nb & Fields2D::SOLID) || (nb & Fields2D::BC_V)) sumP += pij;
        else if (!(nb & Fields2D::AIR)) sumP += d_p[pnx * (j - 1) + i];
    }
    // Top (i, j+1): BC_V checked on current cell
    {
        const uint16_t nb = d_labels[pnx * (j + 1) + i];
        if ((nb & Fields2D::SOLID) || (lbl & Fields2D::BC_V)) sumP += pij;
        else if (!(nb & Fields2D::AIR)) sumP += d_p[pnx * (j + 1) + i];
    }

    const varType rhs = d_rhs[id];

    // True residual before update: r = b + nS - 4*p
    const varType r = rhs + sumP - varType(4) * pij;
    d_resSq[id] = r * r;

    // SOR update
    const varType p_gs = (rhs + sumP) * varType(0.25);
    d_p[id] = pij + omega * (p_gs - pij);
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
    CUDA_CHECK(cudaMalloc(&d.d_cell_idx, nP * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d.d_alive,    nP * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d.d_tmp,      nP * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_perm,     nP * sizeof(int)));

    // ── Grid velocity ──
    CUDA_CHECK(cudaMalloc(&d.d_u,     u_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_v,     v_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_u_num, u_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_u_den, u_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_v_num, v_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_v_den, v_sz * sizeof(varType)));

    // ── Pressure / RBGS workspace ──
    CUDA_CHECK(cudaMalloc(&d.d_p,     p_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_rhs,   p_sz * sizeof(varType)));
    CUDA_CHECK(cudaMalloc(&d.d_resSq, p_sz * sizeof(varType)));
    CUDA_CHECK(cudaMemset(d.d_p,      0, p_sz * sizeof(varType)));

    // ── Labels + cell occupancy ──
    CUDA_CHECK(cudaMalloc(&d.d_labels,     p_sz * sizeof(uint16_t)));
    CUDA_CHECK(cudaMalloc(&d.d_cell_start, c_sz * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d.d_cell_count, c_sz * sizeof(int)));
}

// ─── Destructor ───────────────────────────────────────────────────────────────

GPIC::~GPIC() {
    if (!dev_) return;
    auto& d = *dev_;
#define FREE(p) if (p) { cudaFree(p); (p) = nullptr; }
    FREE(d.d_pos_x)     FREE(d.d_pos_y)
    FREE(d.d_vel_x)     FREE(d.d_vel_y)
    FREE(d.d_cell_idx)  FREE(d.d_alive)
    FREE(d.d_tmp)       FREE(d.d_perm)
    FREE(d.d_u)         FREE(d.d_v)
    FREE(d.d_u_num)     FREE(d.d_u_den)
    FREE(d.d_v_num)     FREE(d.d_v_den)
    FREE(d.d_p)         FREE(d.d_rhs)       FREE(d.d_resSq)
    FREE(d.d_labels)    FREE(d.d_cell_start) FREE(d.d_cell_count)
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
            N, nx, ny, dx, dy,
            d.d_pos_x, d.d_pos_y, d.d_vel_x, d.d_vel_y,
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

// ─── GPU pressure solve ───────────────────────────────────────────────────────
// 1. Compute RHS = -coef*div on GPU
// 2. Run GPU RBGS until ||r||/||b|| < tol
// 3. Apply pressure gradient to d_u, d_v on GPU

void GPIC::GPUMakeIncompressible() {
    auto& d = *dev_;
    const int pnx  = nx + 2;
    const int pny  = ny + 2;
    const int pNxy = pnx * pny;
    const varType coef = d.density * dx * dx / d.dt;

    // ── Compute RHS ──
    CUDA_CHECK(cudaMemset(d.d_rhs,   0, pNxy * sizeof(varType)));
    CUDA_CHECK(cudaMemset(d.d_resSq, 0, pNxy * sizeof(varType)));
    {
        const dim3 block(BX, BY);
        const dim3 grid((nx + BX - 1) / BX, (ny + BY - 1) / BY);
        computeRHSKernel<<<grid, block>>>(
            nx, ny, dx, dy, coef, d.d_labels, d.d_u, d.d_v, d.d_rhs);
        CUDA_CHECK(cudaGetLastError());
    }

    // ── ||b||² for convergence criterion ──
    auto rhs_ptr   = thrust::device_pointer_cast(d.d_rhs);
    const varType bNormSq = thrust::transform_reduce(
        rhs_ptr, rhs_ptr + pNxy, SquareFunctor{}, varType(0), thrust::plus<varType>());
    if (bNormSq < REAL_EPSILON * REAL_EPSILON) return;
    const varType bNorm = static_cast<varType>(std::sqrt(static_cast<double>(bNormSq)));

    // ── RBGS parameters ──
    constexpr double PI = 3.14159265358979323846;
    const int N_min = std::min(nx, ny);
    const varType omega = static_cast<varType>(
        std::min(1.95, 2.0 / (1.0 + std::sin(PI / static_cast<double>(N_min)))));

    const dim3 block2(BX, BY);
    const dim3 grid2((pnx - 2 + BX - 1) / BX, (pny - 2 + BY - 1) / BY);
    const int     maxIters = params.solver.maxIters;
    const varType tol      = static_cast<varType>(params.solver.tolerance);
    auto resSq_ptr = thrust::device_pointer_cast(d.d_resSq);

    // ── RBGS loop (persistent d_p = warm start) ──
    for (int it = 0; it < maxIters; ++it) {
        gpicRBGSKernel<<<grid2, block2>>>(0, pnx, pny, omega,
                                           d.d_labels, d.d_rhs, d.d_p, d.d_resSq);
        gpicRBGSKernel<<<grid2, block2>>>(1, pnx, pny, omega,
                                           d.d_labels, d.d_rhs, d.d_p, d.d_resSq);
        CUDA_CHECK(cudaGetLastError());

        if ((it % CHECK_EVERY) != CHECK_EVERY - 1) continue;

        const varType resNormSq = thrust::reduce(
            resSq_ptr, resSq_ptr + pNxy, varType(0), thrust::plus<varType>());
        if (std::sqrt(static_cast<double>(resNormSq)) <=
            static_cast<double>(tol) * static_cast<double>(bNorm))
            break;
    }

    // ── Apply pressure gradient ──
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
        N, nx, ny, dx, dy, d.dt, params.gravity,
        d.d_u, d.d_v, d.d_pos_x, d.d_pos_y,
        d.d_vel_x, d.d_vel_y);
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

void GPIC::GPUSortParticles() {
    auto& d = *dev_;
    const int N = d.nParticles;

    CUDA_CHECK(cudaMemset(d.d_cell_count, 0,  static_cast<std::size_t>(nx) * ny * sizeof(int)));
    CUDA_CHECK(cudaMemset(d.d_cell_start, -1, static_cast<std::size_t>(nx) * ny * sizeof(int)));

    if (N == 0) return;

    thrust::sequence(thrust::device,
                     thrust::device_pointer_cast(d.d_perm),
                     thrust::device_pointer_cast(d.d_perm + N));

    thrust::sort_by_key(thrust::device,
                        thrust::device_pointer_cast(d.d_cell_idx),
                        thrust::device_pointer_cast(d.d_cell_idx + N),
                        thrust::device_pointer_cast(d.d_perm));

    auto gather_vt = [&](varType* src) {
        thrust::gather(thrust::device,
                       thrust::device_pointer_cast(d.d_perm),
                       thrust::device_pointer_cast(d.d_perm + N),
                       thrust::device_pointer_cast(src),
                       thrust::device_pointer_cast(d.d_tmp));
        CUDA_CHECK(cudaMemcpy(src, d.d_tmp, N * sizeof(varType),
                              cudaMemcpyDeviceToDevice));
    };
    gather_vt(d.d_pos_x);
    gather_vt(d.d_pos_y);
    gather_vt(d.d_vel_x);
    gather_vt(d.d_vel_y);

    const int grid1d = (N + BLOCK - 1) / BLOCK;
    buildCellRangesKernel<<<grid1d, BLOCK>>>(
        N, nx * ny, d.d_cell_idx, d.d_cell_start, d.d_cell_count);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
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
    CUDA_CHECK(cudaDeviceSynchronize());
}

// ─── CFL (GPU) ────────────────────────────────────────────────────────────────

int GPIC::computeAdvectionSubstepsGPU() const {
    auto& d = *dev_;
    const int N = d.nParticles;
    if (N == 0 || params.max_cfl <= REAL_EPSILON || dt <= REAL_EPSILON) return 1;

    auto vx_ptr = thrust::device_pointer_cast(d.d_vel_x);
    auto vy_ptr = thrust::device_pointer_cast(d.d_vel_y);

    const varType max_vx = thrust::transform_reduce(
        vx_ptr, vx_ptr + N, AbsFunctor{}, varType(0), thrust::maximum<varType>());
    const varType max_vy = thrust::transform_reduce(
        vy_ptr, vy_ptr + N, AbsFunctor{}, varType(0), thrust::maximum<varType>());

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
        GPUMakeIncompressible();
        GPUProjectGridOnParticles();
        GPUAdvect();
        GPUSortParticles();
        GPUUpdateCellState();

        if (params.refill) {
            // Refill needs CPU particle inspection
            DownloadAll();
            DownloadParticles();
            const int target = params.ppcx * params.ppcy;
            for (int cj = 0; cj < ny; ++cj) {
                for (int ci = 0; ci < nx; ++ci) {
                    if (!IS_FLUID(fields->Label(ci + 1, cj + 1))) continue;
                    Particles& cell = (*cloud_)(ci, cj);
                    const int missing = target - cell.size();
                    if (missing <= 0) continue;
                    for (int m = 0; m < missing; ++m) {
                        const varType x = (ci + rand01()) * dx;
                        const varType y = (cj + rand01()) * dy;
                        const varType u = fields->u.interpolate<0>(x, y, dx, dy);
                        const varType v = fields->v.interpolate<1>(x, y, dx, dy);
                        cell.Add(x, y, u, v, 0u);
                    }
                }
            }
            dev_->nParticles = cloud_->totalSize();
            UploadParticles();
        }
    }

    setTimeStep(frameDt);
    dev_->dt = frameDt;
}

// ─── Run ──────────────────────────────────────────────────────────────────────

void GPIC::Run() {
    // Upload initial state once — particles and grid stay on GPU throughout
    UploadParticles();
    UploadGrid();

    GPUProjectGridOnParticles();   // init particle velocities from initial grid

    WriteOutput(0);                // downloads from GPU internally
    RunLoop(std::max(1, params.nt / 100));
}

// ─── WriteOutput ──────────────────────────────────────────────────────────────
// Downloads grid from GPU to CPU before calling base-class writers.

void GPIC::WriteOutput(int step) const {
    if (step % params.sampling_rate != 0) return;

    DownloadAll();                          // d_u, d_v, d_p, d_labels → fields
    fields->Div();                          // compute fields->div from u, v
    fields->VelocityNormCenterGrid();       // compute normVelocity from u, v

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
