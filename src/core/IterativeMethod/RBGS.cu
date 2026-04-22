
#ifdef USE_CUDA

#include "IterativeMethods.hpp"
#include "../Fields.hpp"
#include "../Precision.hpp"

#include <cmath>
#include <cstdio>
#include <cuda_runtime.h>

#define CUDA_CHECK(call)                                                        \
    do {                                                                        \
        cudaError_t _e = (call);                                                \
        if (_e != cudaSuccess) {                                                \
            fprintf(stderr, "[CUDA] %s:%d  %s\n",                              \
                    __FILE__, __LINE__, cudaGetErrorString(_e));                \
            std::exit(EXIT_FAILURE);                                            \
        }                                                                       \
    } while (0)

__device__ __forceinline__ int pidx(int pnx, int i, int j) {
    return pnx * j + i;
}
//@todo must be rewritten by ourself more simpler
static constexpr int BX = 32;
static constexpr int BY = 8;

// ── Warp-level max reduction ──────────────────────────────────────────────────
// Butterfly shuffle: each step exchanges values between lanes that differ by
// `mask` bits, keeping the running maximum.  After log2(32)=5 steps every
// lane holds the warp-wide max.
//
// atomicMax has no double overload in CUDA, so we reinterpret the bit pattern
// as uint64.  This is valid for non-negative doubles because IEEE-754 doubles
// with the same sign sort identically as unsigned integers.
__device__ __forceinline__ void warpAtomicMaxRes(double *out, double val) {
    for (int mask = warpSize / 2; mask > 0; mask >>= 1)
        val = fmax(val, __shfl_xor_sync(0xffffffff, val, mask));

    if ((threadIdx.x % warpSize) == 0)
        atomicMax(reinterpret_cast<unsigned long long *>(out),
                  __double_as_longlong(val));
}

// ── RBGS kernel ───────────────────────────────────────────────────────────────
// Each thread owns one interior p-cell.
// A shared-memory tile of (BX+2)×(BY+2) holds centre values + 1-cell halo
// so that all four neighbours are read from fast shared memory.
__global__ void rbgsKernel(
        int    colour,
        int    pnx, int pny,
        double omega,
        const int    * __restrict__ lbl,
        const double * __restrict__ b,      // b[id] = -coef * div(i-1, j-1)
        double       * __restrict__ p,
        double       * __restrict__ max_res_out,
        bool           do_check)
{
    __shared__ double tile[BY + 2][BX + 2];

    const int ti = threadIdx.x;   // 0 .. BX-1
    const int tj = threadIdx.y;   // 0 .. BY-1

    // Map to interior p-space: 1 .. pnx-2, 1 .. pny-2.
    const int gi = blockIdx.x * BX + ti + 1;
    const int gj = blockIdx.y * BY + tj + 1;

    const bool in_domain = (gi <= pnx - 2) && (gj <= pny - 2);

    // ── Load tile (centre + four halos) ──────────────────────────────────────
    tile[tj + 1][ti + 1] = in_domain ? p[pidx(pnx, gi, gj)] : 0.0;

    // Left halo
    if (ti == 0) {
        const int xi = blockIdx.x * BX;          // gi - 1
        tile[tj + 1][0] = (xi >= 1 && xi <= pnx - 2 && gj <= pny - 2)
                          ? p[pidx(pnx, xi, gj)] : 0.0;
    }
    // Right halo
    if (ti == BX - 1) {
        const int xi = blockIdx.x * BX + BX + 1; // gi + 1 of last column
        tile[tj + 1][BX + 1] = (xi >= 1 && xi <= pnx - 2 && gj <= pny - 2)
                                ? p[pidx(pnx, xi, gj)] : 0.0;
    }
    // Bottom halo
    if (tj == 0) {
        const int yj = blockIdx.y * BY;           // gj - 1
        tile[0][ti + 1] = (yj >= 1 && yj <= pny - 2 && gi <= pnx - 2)
                          ? p[pidx(pnx, gi, yj)] : 0.0;
    }
    // Top halo
    if (tj == BY - 1) {
        const int yj = blockIdx.y * BY + BY + 1; // gj + 1 of last row
        tile[BY + 1][ti + 1] = (yj >= 1 && yj <= pny - 2 && gi <= pnx - 2)
                                ? p[pidx(pnx, gi, yj)] : 0.0;
    }
    __syncthreads();

    if (!in_domain)                    return;
    if ((gi + gj) % 2 != colour)      return;

    const int id  = pidx(pnx, gi, gj);
    const int cur = lbl[id];
    if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur)) return;
    // same as the CPP
    const double pij = tile[tj + 1][ti + 1];
    double sumP = 0.0;

    // Left
    { const int nb = lbl[pidx(pnx, gi - 1, gj)];
      if      (IS_SOLID(nb) || IS_BC_U(nb)) sumP += pij;
      else if (!IS_AIR(nb))                 sumP += tile[tj + 1][ti]; }

    // Right
    { const int nb = lbl[pidx(pnx, gi + 1, gj)];
      if      (IS_SOLID(nb) || IS_BC_U(cur)) sumP += pij;
      else if (!IS_AIR(nb))                  sumP += tile[tj + 1][ti + 2]; }

    // Bottom
    { const int nb = lbl[pidx(pnx, gi, gj - 1)];
      if      (IS_SOLID(nb) || IS_BC_V(nb)) sumP += pij;
      else if (!IS_AIR(nb))                 sumP += tile[tj][ti + 1]; }

    // Top
    { const int nb = lbl[pidx(pnx, gi, gj + 1)];
      if      (IS_SOLID(nb) || IS_BC_V(cur)) sumP += pij;
      else if (!IS_AIR(nb))                  sumP += tile[tj + 2][ti + 1]; }
    const double p_gs  = (b[id] + sumP) / 4.0;
    const double p_new = pij + omega * (p_gs - pij);
    p[id] = p_new;

    // Convergence residual
    if (do_check)
        warpAtomicMaxRes(max_res_out, fabs(p_new - pij));
}

bool solveRedBlackGaussSeidel_GPU(Fields2D &fields, int nx, int ny,
                                  double coef, int maxIters, double tol,
                                  double /*beta*/) {
    fields.Div();

    const int pnx = fields.p.nx;   // nx + 2
    const int pny = fields.p.ny;   // ny + 2
    const int N   = pnx * pny;

    // SOR omega: same formula as CPU.
    constexpr double PI = 3.14159265358979323846;
    const int    N_min  = (nx < ny) ? nx : ny;
    const double omega  = std::min(1.95, 2.0 / (1.0 + std::sin(PI / N_min)));

    int    *h_lbl = new int   [N]();
    double *h_b   = new double[N]();
    double *h_p   = new double[N]();
    //@todo to be improved
    for (int j = 0; j < pny; ++j)
        for (int i = 0; i < pnx; ++i)
            h_lbl[pnx * j + i] = static_cast<int>(fields.Label(i, j));

    double b_norm = 0.0;
    for (int j = 1; j < pny - 1; ++j)
        for (int i = 1; i < pnx - 1; ++i) {
            const int id  = pnx * j + i;
            const int cur = h_lbl[id];
            if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur)) continue;
            const double val = -coef * static_cast<double>(
                                   fields.div.Get(i - 1, j - 1));
            h_b[id] = val;
            if (fabs(val) > b_norm) b_norm = fabs(val);
        }

    if (b_norm < 1e-30) {
#ifndef NDEBUG
        printf("  RBGS_GPU: RHS negligible.\n");
#endif
        delete[] h_lbl; delete[] h_b; delete[] h_p;
        return true;
    }

    // Warm-start from current fields.p.
    for (int j = 1; j < pny - 1; ++j)
        for (int i = 1; i < pnx - 1; ++i)
            h_p[pnx * j + i] = static_cast<double>(fields.p.Get(i, j));

    int    *d_lbl; CUDA_CHECK(cudaMalloc(&d_lbl, N * sizeof(int)));
    double *d_b;   CUDA_CHECK(cudaMalloc(&d_b,   N * sizeof(double)));
    double *d_p;   CUDA_CHECK(cudaMalloc(&d_p,   N * sizeof(double)));
    double *d_res; CUDA_CHECK(cudaMalloc(&d_res, sizeof(double)));

    CUDA_CHECK(cudaMemcpy(d_lbl, h_lbl, N * sizeof(int),    cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b,   h_b,   N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_p,   h_p,   N * sizeof(double), cudaMemcpyHostToDevice));

    const dim3 block(BX, BY);
    const dim3 grid((pnx + BX - 1) / BX, (pny + BY - 1) / BY);

    // Pull d_res to host every CHECK_EVERY iters to amortise PCIe latency.
    static constexpr int CHECK_EVERY = 8;
    bool   converged = false;
    double max_res   = 0.0;

    for (int it = 0; it < maxIters; ++it) {
        const bool check = ((it % CHECK_EVERY) == (CHECK_EVERY - 1));

        if (check) {
            const double zero = 0.0;
            CUDA_CHECK(cudaMemcpy(d_res, &zero, sizeof(double),
                                  cudaMemcpyHostToDevice));
        }

        rbgsKernel<<<grid, block>>>(0, pnx, pny, omega,
                                    d_lbl, d_b, d_p, d_res, check);
        rbgsKernel<<<grid, block>>>(1, pnx, pny, omega,
                                    d_lbl, d_b, d_p, d_res, check);

        if (check) {
            CUDA_CHECK(cudaMemcpy(&max_res, d_res, sizeof(double),
                                  cudaMemcpyDeviceToHost));
            if (max_res / b_norm <= tol) {
#ifndef NDEBUG
                printf("  RBGS_GPU converged in %d iters, rel.res = %g\n",
                       it + 1, max_res / b_norm);
#endif
                converged = true;
                break;
            }
        }
    }

    if (!converged)
        printf("  RBGS_GPU: maxIters=%d  rel.res=%g\n",
               maxIters, max_res / b_norm);

    CUDA_CHECK(cudaMemcpy(h_p, d_p, N * sizeof(double), cudaMemcpyDeviceToHost));

    for (int j = 1; j < pny - 1; ++j)
        for (int i = 1; i < pnx - 1; ++i) {
            const int cur = h_lbl[pnx * j + i];
            if (!IS_SOLID(cur) && !IS_AIR(cur) && !IS_BC_P(cur))
                fields.p.Set(i, j, static_cast<varType>(h_p[pnx * j + i]));
        }

    CUDA_CHECK(cudaFree(d_lbl));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_p));
    CUDA_CHECK(cudaFree(d_res));
    delete[] h_lbl; delete[] h_b; delete[] h_p;

    return converged;
}

#endif // USE_CUDA
