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
static constexpr int CHECK_EVERY = 32;

__host__ __device__ __forceinline__ int pidx(int pnx, int i, int j) {
  return pnx * j + i;
}

// Kernel uses varType so it works at native precision (float or double).
__global__ void rbgsColorKernel(int colour, int pnx, int pny, varType omega,
                                const int *lbl, const varType *b, varType *p,
                                varType *deltaSq, bool accumulateResidual) {
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

  const varType pij = p[id];
  varType sumP = varType(0);

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

  const varType p_gs = (b[id] + sumP) * varType(0.25);
  p[id] = pij + omega * (p_gs - pij);

  if (accumulateResidual) {
    // True residual before update: r = b + nS - 4*p
    const varType r = b[id] + sumP - varType(4) * pij;
    deltaSq[id] = r * r;
  }
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
  // Compute optimal SOR omega in double for accuracy, then cast.
  const varType omega = static_cast<varType>(
      std::min(1.95, 2.0 / (1.0 + std::sin(PI / static_cast<double>(N_min)))));

  std::vector<int>     h_lbl(static_cast<std::size_t>(N), 0);
  std::vector<varType> h_b(static_cast<std::size_t>(N), varType(0));
  std::vector<varType> h_p(static_cast<std::size_t>(N), varType(0));
  std::vector<varType> h_deltaSq(static_cast<std::size_t>(N), varType(0));

  for (int j = 0; j < pny; ++j) {
    for (int i = 0; i < pnx; ++i) {
      h_lbl[pidx(pnx, i, j)] = static_cast<int>(fields.Label(i, j));
      h_p[pidx(pnx, i, j)]   = fields.p.Get(i, j);
    }
  }

  // Accumulate bNormSq in double to avoid catastrophic cancellation.
  double bNormSq = 0.0;
  for (int j = 1; j < pny - 1; ++j) {
    for (int i = 1; i < pnx - 1; ++i) {
      const int id  = pidx(pnx, i, j);
      const int cur = h_lbl[id];
      if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur))
        continue;

      const varType b_val =
          -coef * fields.div.Get(i - 1, j - 1);
      h_b[id] = b_val;
      bNormSq += static_cast<double>(b_val) * static_cast<double>(b_val);
    }
  }

  if (bNormSq < 1e-60) {
    DBG_PRINTF("RBGS_GPU: no fluid pressure cells or zero RHS");
    return true;
  }
  const double bNorm = std::sqrt(bNormSq);

  int     *d_lbl     = nullptr;
  varType *d_b       = nullptr;
  varType *d_p       = nullptr;
  varType *d_deltaSq = nullptr;

  CUDA_CHECK(cudaMalloc(&d_lbl,     static_cast<std::size_t>(N) * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&d_b,       static_cast<std::size_t>(N) * sizeof(varType)));
  CUDA_CHECK(cudaMalloc(&d_p,       static_cast<std::size_t>(N) * sizeof(varType)));
  CUDA_CHECK(cudaMalloc(&d_deltaSq, static_cast<std::size_t>(N) * sizeof(varType)));

  CUDA_CHECK(cudaMemcpy(d_lbl, h_lbl.data(),
                        static_cast<std::size_t>(N) * sizeof(int),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_b, h_b.data(),
                        static_cast<std::size_t>(N) * sizeof(varType),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_p, h_p.data(),
                        static_cast<std::size_t>(N) * sizeof(varType),
                        cudaMemcpyHostToDevice));

  const dim3 block(BX, BY);
  const dim3 grid((pnx - 2 + BX - 1) / BX, (pny - 2 + BY - 1) / BY);

  bool converged = false;

  for (int it = 0; it < maxIters; ++it) {
    const bool check = ((it % CHECK_EVERY) == (CHECK_EVERY - 1));
    if (check) {
      CUDA_CHECK(cudaMemset(d_deltaSq, 0,
                            static_cast<std::size_t>(N) * sizeof(varType)));
    }

    rbgsColorKernel<<<grid, block>>>(0, pnx, pny, omega, d_lbl, d_b, d_p,
                                     d_deltaSq, check);
    rbgsColorKernel<<<grid, block>>>(1, pnx, pny, omega, d_lbl, d_b, d_p,
                                     d_deltaSq, check);
    CUDA_CHECK(cudaGetLastError());

    if (!check)
      continue;

    CUDA_CHECK(cudaMemcpy(h_deltaSq.data(), d_deltaSq,
                          static_cast<std::size_t>(N) * sizeof(varType),
                          cudaMemcpyDeviceToHost));

    // Accumulate in double for accurate convergence check.
    double sumSq = 0.0;
    for (varType v : h_deltaSq)
      sumSq += static_cast<double>(v);

    const double relRes = std::sqrt(sumSq) / bNorm;
    if (relRes <= static_cast<double>(tol)) {
      DBG_PRINTF("RBGS_GPU converged in %d iters, rel.res = %.6g", it + 1, relRes);
      converged = true;
      break;
    }
  }

  if (!converged)
    DBG_PRINTF("RBGS_GPU: reached maxIters = %d", maxIters);

  CUDA_CHECK(cudaMemcpy(h_p.data(), d_p,
                        static_cast<std::size_t>(N) * sizeof(varType),
                        cudaMemcpyDeviceToHost));

  for (int j = 1; j < pny - 1; ++j) {
    for (int i = 1; i < pnx - 1; ++i) {
      const int id  = pidx(pnx, i, j);
      const int cur = h_lbl[id];
      if (!IS_SOLID(cur) && !IS_AIR(cur) && !IS_BC_P(cur))
        fields.p.Set(i, j, h_p[id]);
    }
  }

  CUDA_CHECK(cudaFree(d_lbl));
  CUDA_CHECK(cudaFree(d_b));
  CUDA_CHECK(cudaFree(d_p));
  CUDA_CHECK(cudaFree(d_deltaSq));

  return converged;
}

#endif // USE_CUDA
