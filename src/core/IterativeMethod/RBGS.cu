#ifdef USE_CUDA

#include "IterativeMethods.hpp"
#include "../Fields.hpp"
#include "../Precision.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <vector>

#define CUDA_CHECK(call)                                                        \
  do {                                                                          \
    cudaError_t _e = (call);                                                    \
    if (_e != cudaSuccess) {                                                    \
      fprintf(stderr, "[CUDA] %s:%d  %s\n", __FILE__, __LINE__,                 \
              cudaGetErrorString(_e));                                          \
      std::exit(EXIT_FAILURE);                                                  \
    }                                                                           \
  } while (0)

namespace {

static constexpr int BX = 32;
static constexpr int BY = 8;
static constexpr int CHECK_EVERY = 8;

__host__ __device__ __forceinline__ int pidx(int pnx, int i, int j) {
  return pnx * j + i;
}

__global__ void rbgsColorKernel(int colour, int pnx, int pny, double omega,
                                const int *lbl, const double *b, double *p,
                                double *deltaSq, bool accumulateResidual) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  const int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

  if (i >= pnx - 1 || j >= pny - 1)
    return;
  if (((i + j) & 1) != colour)
    return;

  const int id = pidx(pnx, i, j);
  const int cur = lbl[id];
  if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur))
    return;

  const double pij = p[id];
  double sumP = 0.0;

  // Left.
  {
    const int nb = lbl[pidx(pnx, i - 1, j)];
    if (IS_SOLID(nb) || IS_BC_U(nb))
      sumP += pij;
    else if (!IS_AIR(nb))
      sumP += p[pidx(pnx, i - 1, j)];
  }

  // Right.
  {
    const int nb = lbl[pidx(pnx, i + 1, j)];
    if (IS_SOLID(nb) || IS_BC_U(cur))
      sumP += pij;
    else if (!IS_AIR(nb))
      sumP += p[pidx(pnx, i + 1, j)];
  }

  // Bottom.
  {
    const int nb = lbl[pidx(pnx, i, j - 1)];
    if (IS_SOLID(nb) || IS_BC_V(nb))
      sumP += pij;
    else if (!IS_AIR(nb))
      sumP += p[pidx(pnx, i, j - 1)];
  }

  // Top.
  {
    const int nb = lbl[pidx(pnx, i, j + 1)];
    if (IS_SOLID(nb) || IS_BC_V(cur))
      sumP += pij;
    else if (!IS_AIR(nb))
      sumP += p[pidx(pnx, i, j + 1)];
  }

  const double p_gs = (b[id] + sumP) / 4.0;
  const double p_new = pij + omega * (p_gs - pij);
  p[id] = p_new;

  if (accumulateResidual)
    deltaSq[id] = (p_new - pij) * (p_new - pij);
}

} // namespace

bool solveRedBlackGaussSeidel_GPU(Fields2D &fields, int nx, int ny,
                                  varType coef, int maxIters, varType tol,
                                  varType /*beta*/) {
  fields.Div();

  const int pnx = fields.p.nx;
  const int pny = fields.p.ny;
  const int N = pnx * pny;

  constexpr double PI = 3.14159265358979323846;
  const int N_min = (nx < ny) ? nx : ny;
  const double omega =
      std::min(1.95, 2.0 / (1.0 + std::sin(PI / static_cast<double>(N_min))));

  std::vector<int> h_lbl(static_cast<std::size_t>(N), 0);
  std::vector<double> h_b(static_cast<std::size_t>(N), 0.0);
  std::vector<double> h_p(static_cast<std::size_t>(N), 0.0);
  std::vector<double> h_deltaSq(static_cast<std::size_t>(N), 0.0);

  int fluidCount = 0;
  for (int j = 0; j < pny; ++j) {
    for (int i = 0; i < pnx; ++i) {
      h_lbl[pidx(pnx, i, j)] = static_cast<int>(fields.Label(i, j));
      h_p[pidx(pnx, i, j)] = static_cast<double>(fields.p.Get(i, j));
    }
  }

  for (int j = 1; j < pny - 1; ++j) {
    for (int i = 1; i < pnx - 1; ++i) {
      const int id = pidx(pnx, i, j);
      const int cur = h_lbl[id];
      if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur))
        continue;

      h_b[id] = -coef * static_cast<double>(fields.div.Get(i - 1, j - 1));
      ++fluidCount;
    }
  }

  if (fluidCount == 0) {
    DBG_PRINTF("RBGS_GPU: no fluid pressure cells");
    return true;
  }

  int *d_lbl = nullptr;
  double *d_b = nullptr;
  double *d_p = nullptr;
  double *d_deltaSq = nullptr;

  CUDA_CHECK(cudaMalloc(&d_lbl, static_cast<std::size_t>(N) * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&d_b, static_cast<std::size_t>(N) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_p, static_cast<std::size_t>(N) * sizeof(double)));
  CUDA_CHECK(
      cudaMalloc(&d_deltaSq, static_cast<std::size_t>(N) * sizeof(double)));

  CUDA_CHECK(cudaMemcpy(d_lbl, h_lbl.data(), static_cast<std::size_t>(N) *
                                            sizeof(int),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_b, h_b.data(), static_cast<std::size_t>(N) *
                                          sizeof(double),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_p, h_p.data(), static_cast<std::size_t>(N) *
                                          sizeof(double),
                        cudaMemcpyHostToDevice));

  const dim3 block(BX, BY);
  const dim3 grid((pnx - 2 + BX - 1) / BX, (pny - 2 + BY - 1) / BY);

  bool converged = false;
  double res = 0.0;
  double res0 = -1.0;

  for (int it = 0; it < maxIters; ++it) {
    const bool check = ((it % CHECK_EVERY) == (CHECK_EVERY - 1));
    if (check) {
      CUDA_CHECK(cudaMemset(d_deltaSq, 0,
                            static_cast<std::size_t>(N) * sizeof(double)));
    }

    rbgsColorKernel<<<grid, block>>>(0, pnx, pny, omega, d_lbl, d_b, d_p,
                                     d_deltaSq, check);
    rbgsColorKernel<<<grid, block>>>(1, pnx, pny, omega, d_lbl, d_b, d_p,
                                     d_deltaSq, check);
    CUDA_CHECK(cudaGetLastError());

    if (!check)
      continue;

    CUDA_CHECK(cudaMemcpy(h_deltaSq.data(), d_deltaSq,
                          static_cast<std::size_t>(N) * sizeof(double),
                          cudaMemcpyDeviceToHost));

    double sumSq = 0.0;
    for (double v : h_deltaSq)
      sumSq += v;

    res = std::sqrt(sumSq / static_cast<double>(fluidCount));
    if (res0 < 0.0) {
      res0 = res;
      if (res0 < 1e-30) {
        converged = true;
        break;
      }
      continue;
    }

    if (res0 < 1e-30 || res / res0 <= tol) {
      DBG_PRINTF("RBGS_GPU converged in %d iters, rel.res = %.6g", it + 1,
                 res / res0);
      converged = true;
      break;
    }
  }

  if (!converged)
    DBG_PRINTF("RBGS_GPU: reached maxIters = %d", maxIters);

  CUDA_CHECK(cudaMemcpy(h_p.data(), d_p, static_cast<std::size_t>(N) *
                                          sizeof(double),
                        cudaMemcpyDeviceToHost));

  for (int j = 1; j < pny - 1; ++j) {
    for (int i = 1; i < pnx - 1; ++i) {
      const int id = pidx(pnx, i, j);
      const int cur = h_lbl[id];
      if (!IS_SOLID(cur) && !IS_AIR(cur) && !IS_BC_P(cur))
        fields.p.Set(i, j, static_cast<varType>(h_p[id]));
    }
  }

  CUDA_CHECK(cudaFree(d_lbl));
  CUDA_CHECK(cudaFree(d_b));
  CUDA_CHECK(cudaFree(d_p));
  CUDA_CHECK(cudaFree(d_deltaSq));

  return converged;
}

#endif // USE_CUDA
