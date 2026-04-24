#ifdef USE_CUDA

#include "IterativeMethods.hpp"
#include "../Fields.hpp"
#include "../Precision.hpp"

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

static constexpr int THREADS = 256;

__host__ __device__ __forceinline__ int pidx(int pnx, int i, int j) {
  return pnx * j + i;
}

__global__ void applyAKernel(int N, int pnx, int pny, const int *lbl,
                             const double *x, double *Ax) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id >= N)
    return;

  const int i = id % pnx;
  const int j = id / pnx;
  if (i == 0 || i == pnx - 1 || j == 0 || j == pny - 1) {
    Ax[id] = 0.0;
    return;
  }

  const int cur = lbl[id];
  if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur)) {
    Ax[id] = 0.0;
    return;
  }

  double acc = 4.0 * x[id];

  // Left.
  {
    const int nb = lbl[pidx(pnx, i - 1, j)];
    if (IS_SOLID(nb) || IS_BC_U(nb))
      acc -= x[id];
    else if (!IS_AIR(nb) && !IS_BC_P(nb))
      acc -= x[pidx(pnx, i - 1, j)];
  }

  // Right.
  {
    const int nb = lbl[pidx(pnx, i + 1, j)];
    if (IS_SOLID(nb) || IS_BC_U(cur))
      acc -= x[id];
    else if (!IS_AIR(nb) && !IS_BC_P(nb))
      acc -= x[pidx(pnx, i + 1, j)];
  }

  // Bottom.
  {
    const int nb = lbl[pidx(pnx, i, j - 1)];
    if (IS_SOLID(nb) || IS_BC_V(nb))
      acc -= x[id];
    else if (!IS_AIR(nb) && !IS_BC_P(nb))
      acc -= x[pidx(pnx, i, j - 1)];
  }

  // Top.
  {
    const int nb = lbl[pidx(pnx, i, j + 1)];
    if (IS_SOLID(nb) || IS_BC_V(cur))
      acc -= x[id];
    else if (!IS_AIR(nb) && !IS_BC_P(nb))
      acc -= x[pidx(pnx, i, j + 1)];
  }

  Ax[id] = acc;
}

__global__ void initResidualKernel(int N, const int *lbl, const double *b,
                                   const double *Ax, double *r, double *d) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id >= N)
    return;

  const int cur = lbl[id];
  if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur)) {
    r[id] = 0.0;
    d[id] = 0.0;
    return;
  }

  const double val = b[id] - Ax[id];
  r[id] = val;
  d[id] = val;
}

__global__ void updateXandRKernel(int N, const int *lbl, double alpha,
                                  const double *d, const double *q, double *x,
                                  double *r) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id >= N)
    return;

  const int cur = lbl[id];
  if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur))
    return;

  x[id] += alpha * d[id];
  r[id] -= alpha * q[id];
}

__global__ void updateSearchKernel(int N, const int *lbl, double beta,
                                   const double *r, double *d) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id >= N)
    return;

  const int cur = lbl[id];
  if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur)) {
    d[id] = 0.0;
    return;
  }

  d[id] = r[id] + beta * d[id];
}

__global__ void dotKernel(const double *a, const double *b, double *partial,
                          int N) {
  __shared__ double cache[THREADS];

  const int tid = threadIdx.x;
  int id = blockIdx.x * blockDim.x + tid;
  const int stride = blockDim.x * gridDim.x;

  double sum = 0.0;
  while (id < N) {
    sum += a[id] * b[id];
    id += stride;
  }

  cache[tid] = sum;
  __syncthreads();

  for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
    if (tid < offset)
      cache[tid] += cache[tid + offset];
    __syncthreads();
  }

  if (tid == 0)
    partial[blockIdx.x] = cache[0];
}

double dotDevice(const double *a, const double *b, double *d_partial,
                 std::vector<double> &h_partial, int N) {
  if (N <= 0)
    return 0.0;

  const int blocks = static_cast<int>(h_partial.size());
  dotKernel<<<blocks, THREADS>>>(a, b, d_partial, N);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaMemcpy(h_partial.data(), d_partial,
                        static_cast<std::size_t>(blocks) * sizeof(double),
                        cudaMemcpyDeviceToHost));

  double sum = 0.0;
  for (double v : h_partial)
    sum += v;
  return sum;
}

} // namespace

bool solveCG_GPU(Fields2D &fields, varType coef, varType /*beta*/,
                 int maxIters, varType tol) {
  fields.Div();

  const int pnx = fields.p.nx;
  const int pny = fields.p.ny;
  const int N = pnx * pny;

  std::vector<int> h_lbl(static_cast<std::size_t>(N), 0);
  std::vector<double> h_b(static_cast<std::size_t>(N), 0.0);
  std::vector<double> h_x(static_cast<std::size_t>(N), 0.0);

  double bNorm2Sq = 0.0;
  for (int j = 0; j < pny; ++j) {
    for (int i = 0; i < pnx; ++i) {
      const int id = pidx(pnx, i, j);
      h_lbl[id] = static_cast<int>(fields.Label(i, j));
      h_x[id] = static_cast<double>(fields.p.Get(i, j));
    }
  }

  for (int j = 1; j < pny - 1; ++j) {
    for (int i = 1; i < pnx - 1; ++i) {
      const int id = pidx(pnx, i, j);
      const int cur = h_lbl[id];
      if (IS_SOLID(cur) || IS_AIR(cur) || IS_BC_P(cur))
        continue;

      double rhs = -static_cast<double>(coef) *
                   static_cast<double>(fields.div.Get(i - 1, j - 1));

      // Known pressure neighbours contribute directly to the RHS.
      {
        const int nb = h_lbl[pidx(pnx, i - 1, j)];
        if (!(IS_SOLID(nb) || IS_BC_U(nb)) && IS_BC_P(nb))
          rhs += static_cast<double>(fields.p.Get(i - 1, j));
      }
      {
        const int nb = h_lbl[pidx(pnx, i + 1, j)];
        if (!(IS_SOLID(nb) || IS_BC_U(cur)) && IS_BC_P(nb))
          rhs += static_cast<double>(fields.p.Get(i + 1, j));
      }
      {
        const int nb = h_lbl[pidx(pnx, i, j - 1)];
        if (!(IS_SOLID(nb) || IS_BC_V(nb)) && IS_BC_P(nb))
          rhs += static_cast<double>(fields.p.Get(i, j - 1));
      }
      {
        const int nb = h_lbl[pidx(pnx, i, j + 1)];
        if (!(IS_SOLID(nb) || IS_BC_V(cur)) && IS_BC_P(nb))
          rhs += static_cast<double>(fields.p.Get(i, j + 1));
      }

      h_b[id] = rhs;
      bNorm2Sq += rhs * rhs;
    }
  }

  const double bNorm2 = std::sqrt(bNorm2Sq);
  if (bNorm2 < 1e-30) {
    DBG_PRINTF("CG_GPU: RHS negligible");
    return true;
  }

  const int blocks = std::max(1, (N + THREADS - 1) / THREADS);
  std::vector<double> h_partial(static_cast<std::size_t>(blocks), 0.0);

  int *d_lbl = nullptr;
  double *d_b = nullptr;
  double *d_x = nullptr;
  double *d_r = nullptr;
  double *d_d = nullptr;
  double *d_q = nullptr;
  double *d_partial = nullptr;

  CUDA_CHECK(cudaMalloc(&d_lbl, static_cast<std::size_t>(N) * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&d_b, static_cast<std::size_t>(N) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x, static_cast<std::size_t>(N) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_r, static_cast<std::size_t>(N) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_d, static_cast<std::size_t>(N) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_q, static_cast<std::size_t>(N) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_partial,
                        static_cast<std::size_t>(blocks) * sizeof(double)));

  CUDA_CHECK(cudaMemcpy(d_lbl, h_lbl.data(), static_cast<std::size_t>(N) *
                                            sizeof(int),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_b, h_b.data(), static_cast<std::size_t>(N) *
                                          sizeof(double),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_x, h_x.data(), static_cast<std::size_t>(N) *
                                          sizeof(double),
                        cudaMemcpyHostToDevice));

  applyAKernel<<<blocks, THREADS>>>(N, pnx, pny, d_lbl, d_x, d_q);
  CUDA_CHECK(cudaGetLastError());
  initResidualKernel<<<blocks, THREADS>>>(N, d_lbl, d_b, d_q, d_r, d_d);
  CUDA_CHECK(cudaGetLastError());

  double sigma = dotDevice(d_r, d_r, d_partial, h_partial, N);
  if (sigma < 1e-30) {
    DBG_PRINTF("CG_GPU: initial residual negligible");
    CUDA_CHECK(cudaMemcpy(h_x.data(), d_x, static_cast<std::size_t>(N) *
                                            sizeof(double),
                          cudaMemcpyDeviceToHost));
    for (int j = 1; j < pny - 1; ++j)
      for (int i = 1; i < pnx - 1; ++i) {
        const int id = pidx(pnx, i, j);
        const int cur = h_lbl[id];
        if (!IS_SOLID(cur) && !IS_AIR(cur) && !IS_BC_P(cur))
          fields.p.Set(i, j, static_cast<varType>(h_x[id]));
      }

    CUDA_CHECK(cudaFree(d_lbl));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_r));
    CUDA_CHECK(cudaFree(d_d));
    CUDA_CHECK(cudaFree(d_q));
    CUDA_CHECK(cudaFree(d_partial));
    return true;
  }

  bool converged = false;

  for (int it = 0; it < maxIters; ++it) {
    applyAKernel<<<blocks, THREADS>>>(N, pnx, pny, d_lbl, d_d, d_q);
    CUDA_CHECK(cudaGetLastError());

    const double dAd = dotDevice(d_d, d_q, d_partial, h_partial, N);
    if (!std::isfinite(dAd) || std::abs(dAd) < 1e-30)
      break;

    const double alpha = sigma / dAd;
    updateXandRKernel<<<blocks, THREADS>>>(N, d_lbl, alpha, d_d, d_q, d_x, d_r);
    CUDA_CHECK(cudaGetLastError());

    const double sigmaNew = dotDevice(d_r, d_r, d_partial, h_partial, N);
    const double relRes = std::sqrt(sigmaNew) / bNorm2;
    if (relRes <= static_cast<double>(tol)) {
      DBG_PRINTF("CG_GPU converged in %d iters, rel.res = %.6g", it + 1,
                 relRes);
      converged = true;
      break;
    }

    if (!std::isfinite(sigmaNew) || sigmaNew < 1e-30)
      break;

    const double betaCG = sigmaNew / sigma;
    updateSearchKernel<<<blocks, THREADS>>>(N, d_lbl, betaCG, d_r, d_d);
    CUDA_CHECK(cudaGetLastError());
    sigma = sigmaNew;
  }

  if (!converged)
    DBG_PRINTF("CG_GPU: reached maxIters = %d", maxIters);

  CUDA_CHECK(cudaMemcpy(h_x.data(), d_x, static_cast<std::size_t>(N) *
                                          sizeof(double),
                        cudaMemcpyDeviceToHost));
  for (int j = 1; j < pny - 1; ++j) {
    for (int i = 1; i < pnx - 1; ++i) {
      const int id = pidx(pnx, i, j);
      const int cur = h_lbl[id];
      if (!IS_SOLID(cur) && !IS_AIR(cur) && !IS_BC_P(cur))
        fields.p.Set(i, j, static_cast<varType>(h_x[id]));
    }
  }

  CUDA_CHECK(cudaFree(d_lbl));
  CUDA_CHECK(cudaFree(d_b));
  CUDA_CHECK(cudaFree(d_x));
  CUDA_CHECK(cudaFree(d_r));
  CUDA_CHECK(cudaFree(d_d));
  CUDA_CHECK(cudaFree(d_q));
  CUDA_CHECK(cudaFree(d_partial));

  return converged;
}

#endif // USE_CUDA
