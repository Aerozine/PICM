#include "IterativeMethods.hpp"
#include "Fields.hpp"
#include "Grid2D.hpp"
#include "Precision.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
  // sumP,nb

  /*[[nodiscard]] inline std::pair<varType, int>  neighbourSum(const Fields2D &f,
                                              int nx, int ny, int i, int j, double beta) {
  // all 4 neighbours nearly always exist.
  // No bounds checks, no branches.
  if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1) {
    // early exit
    return {f.p.Get(i + 1, j) + f.p.Get(i - 1, j) + f.p.Get(i, j + 1) +
                f.p.Get(i, j - 1),
            4};
  }
  // Border cells only
  double sumP = 0.0;
  int nb = 0;
  if (i + 1 < nx) {sumP += f.p.Get(i + 1, j);++nb;}
  if (i - 1 >= 0) {sumP += f.p.Get(i - 1, j);++nb;}
  if (j + 1 < ny) {sumP += f.p.Get(i, j + 1);++nb;}
  if (j - 1 >= 0) {sumP += f.p.Get(i, j - 1);++nb;}
  return {sumP, nb};
}*/

  [[nodiscard]] inline std::pair<varType, int> neighbourSum(const Fields2D &f,
                     int nx, int ny, int i, int j, double beta) {
  
  if (f.Label(i, j) & Fields2D::SOLID) {
    return {f.p.Get(i, j), -1};
  }

  varType sumP = 0.0;
  int nb = 0;
  const varType pC = f.p.Get(i, j);

  // Left neighbour
  if (i - 1 >= 0) {
    if (f.Label(i - 1, j) & Fields2D::SOLID) {
      sumP += pC - beta * f.u.Get(i, j);
    } else {
      sumP += f.p.Get(i - 1, j);
    }
    ++nb;
  }

  // Right neighbour
  if (i + 1 < nx) {
    if (f.Label(i + 1, j) & Fields2D::SOLID) {
      sumP += pC + beta * f.u.Get(i + 1, j);
    } else {
      sumP += f.p.Get(i + 1, j);
    }
    ++nb;
  }

  // Bottom neighbour
  if (j - 1 >= 0) {
    if (f.Label(i, j - 1) & Fields2D::SOLID) {
      sumP += pC - beta * f.v.Get(i, j);
    } else {
      sumP += f.p.Get(i, j - 1);
    }
    ++nb;
  }

  // Top neighbour
  if (j + 1 < ny) {
    if (f.Label(i, j + 1) & Fields2D::SOLID) {
      sumP += pC + beta * f.v.Get(i, j + 1);
    } else {
      sumP += f.p.Get(i, j + 1);
    }
    ++nb;
  }

  return {sumP, nb};
}

// Gauss-Seidel update for a single FLUID cell.
//  p_new = ( -coef * div_{ij} + sum p_nb ) / N_nb
[[nodiscard]] inline double gsUpdate(const Fields2D &f, int nx, int ny, int i,
                                     int j, double coef, double beta) {
  if (f.Label(i, j) & Fields2D::SOLID)
    return f.p.Get(i,j);
  const auto [sumP, nb] = neighbourSum(f, nx, ny, i, j, beta);
  // if there is no neighbour just keep the same value
  return nb > 0 ? (-coef * f.div.Get(i, j) + sumP) / nb : f.p.Get(i,j);
}

// RMS residual  r_{ij} = -coef·div_{ij} − (N_nb·p_{ij} − Σ p_nb)
// over all FLUID cells.
[[nodiscard]] double residualNorm(const Fields2D &f, int nx, int ny,
                                  double coef, double beta) {
  double sumSq = 0.0;
  int count = 0;

OMP_PRAGMA(omp parallel for collapse(2) reduction(+:sumSq) reduction(+:count))
for (int j = 0; j < ny; ++j) {
  for (int i = 0; i < nx; ++i) {
    if (f.Label(i, j) != Fields2D::FLUID)
      continue;
    const auto [sumP, nb] = neighbourSum(f, nx, ny, i, j, beta);
    const double r = (-coef * f.div.Get(i, j)) - (nb * f.p.Get(i, j) - sumP);
    sumSq += r * r;
    ++count;
  }
}
return (count > 0) ? std::sqrt(sumSq / count) : 0.0;
}

// Relative convergence: records res0 on first call (it==0).
inline bool checkConvergence(double res, double &res0, int it, double tol) {
  if (it == 0) {
    res0 = res;
    return res0 < 1e-30;
  }
  return (res / res0) < tol;
}

// For MICCG0 , multiply the coefficient matrix A times a vector
inline void applyA(const Fields2D &f, int nx, int ny, double scale,
                   const std::vector<double> &Adiag,
                   const std::vector<double> &x, std::vector<double> &y) {

OMP_PRAGMA(omp parallel for collapse(2))
for (int j = 0; j < ny; ++j) {
  for (int i = 0; i < nx; ++i) {
    const int id = nx * j + i;
    if (f.Label(i, j) != Fields2D::FLUID) {
      y[id] = 0.0;
      continue;
    }

    double val = Adiag[id] * x[id];
    if (i + 1 < nx && f.Label(i + 1, j) == Fields2D::FLUID)
      val -= scale * x[nx * j + (i + 1)];
    if (j + 1 < ny && f.Label(i, j + 1) == Fields2D::FLUID)
      val -= scale * x[nx * (j + 1) + i];
    if (i - 1 >= 0 && f.Label(i - 1, j) == Fields2D::FLUID)
      val -= scale * x[nx * j + (i - 1)];
    if (j - 1 >= 0 && f.Label(i, j - 1) == Fields2D::FLUID)
      val -= scale * x[nx * (j - 1) + i];
    y[id] = val;
  }
}
}

// Build precon[id] = 1/E[id]
// CF ugly formula at 5.7
void buildPrecon(const Fields2D &f, int nx, int ny, double scale,
                 const std::vector<double> &Adiag,
                 std::vector<double> &precon) {
  // recommanded value
  constexpr double tau = 0.97;
  constexpr double sigma = 0.25;
  const double scale2 = scale * scale;
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (f.Label(i, j) != Fields2D::FLUID)
        continue;
      const int id = nx * j + i;
      double e = Adiag[id];

      // x lower neighbour (i-1, j)
      if (i - 1 >= 0 && f.Label(i - 1, j) == Fields2D::FLUID) {
        const double pre = precon[nx * j + (i - 1)];
        const double pre2 = pre * pre;
        e -= scale2 * pre2;
        if (j + 1 < ny && f.Label(i - 1, j + 1) == Fields2D::FLUID)
          e -= tau * scale2 * pre2;
      }

      // y lower neighbour (i, j-1)
      if (j - 1 >= 0 && f.Label(i, j - 1) == Fields2D::FLUID) {
        const double pre = precon[nx * (j - 1) + i];
        const double pre2 = pre * pre;
        e -= scale2 * pre2;
        if (i + 1 < nx && f.Label(i + 1, j - 1) == Fields2D::FLUID)
          e -= tau * scale2 * pre2;
      }

      if (e < sigma * Adiag[id])
        e = Adiag[id];
      precon[id] = 1.0 / std::sqrt(e);
    }
  }
}

// CF figure 5.8
void applyPrecon(const Fields2D &f, int nx, int ny, double scale,
                 const std::vector<double> &precon,
                 const std::vector<double> &r, std::vector<double> &z,
                 std::vector<double> &q) {
  // Lq=b forward substitution
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (f.Label(i, j) != Fields2D::FLUID)
        continue;
      const int id = nx * j + i;
      double t = r[id];
      if (i - 1 >= 0 && f.Label(i - 1, j) == Fields2D::FLUID)
        t += scale * precon[nx * j + (i - 1)] * q[nx * j + (i - 1)];
      if (j - 1 >= 0 && f.Label(i, j - 1) == Fields2D::FLUID)
        t += scale * precon[nx * (j - 1) + i] * q[nx * (j - 1) + i];
      q[id] = t * precon[id];
    }
  }

  // L^T p = q backward substitution
  std::fill(z.begin(), z.end(), 0.0);
  for (int j = ny - 1; j >= 0; --j) {
    for (int i = nx - 1; i >= 0; --i) {
      if (f.Label(i, j) != Fields2D::FLUID)
        continue;
      const int id = nx * j + i;
      double t = q[id];
      if (i + 1 < nx && f.Label(i + 1, j) == Fields2D::FLUID)
        t += scale * precon[id] * z[nx * j + (i + 1)];
      if (j + 1 < ny && f.Label(i, j + 1) == Fields2D::FLUID)
        t += scale * precon[id] * z[nx * (j + 1) + i];
      z[id] = t * precon[id];
    }
  }
}

}

void solveJacobi(Fields2D &fields, int nx, int ny, double coef, int maxIters,
                 double tol, double beta) {
  fields.Div();
  Grid2D pNew(nx, ny);
  double res0 = 1.0;

  for (int it = 0; it < maxIters; ++it) {
OMP_PRAGMA(omp parallel for collapse(2))
for (int j = 0; j < ny; ++j)
  for (int i = 0; i < nx; ++i)
    pNew.Set(i, j, gsUpdate(fields, nx, ny, i, j, coef, beta));

OMP_PRAGMA(omp parallel for collapse(2))
for (int j = 0; j < ny; ++j)
  for (int i = 0; i < nx; ++i)
    if (fields.Label(i, j) == Fields2D::FLUID)
      fields.p.Set(i, j, pNew.Get(i, j));

const double res = residualNorm(fields, nx, ny, coef, beta);
if (checkConvergence(res, res0, it, tol)) {
#ifndef NDEBUG
  std::cout << "  Jacobi converged in " << it + 1
            << " iters, rel.res = " << res / res0 << '\n';
#endif
  return;
}
  }
#ifndef NDEBUG
  std::cout << "  Jacobi: reached maxIters = " << maxIters << '\n';
#endif
}

void solveGaussSeidel(Fields2D &fields, int nx, int ny, double coef,
                      int maxIters, double tol, double beta) {
  fields.Div();
  double res0 = 1.0;

  for (int it = 0; it < maxIters; ++it) {
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
          fields.p.Set(i, j, gsUpdate(fields,nx,ny,i,j,coef, beta));
      }

    const double res = residualNorm(fields, nx, ny, coef, beta);
    if (checkConvergence(res, res0, it, tol)) {
#ifndef NDEBUG
      std::cout << "  GaussSeidel converged in " << it + 1
                << " iters, rel.res = " << res / res0 << '\n';
#endif
      return;
    }
  }
#ifndef NDEBUG
  std::cout << "  GaussSeidel: reached maxIters = " << maxIters << '\n';
#endif
}

void solveRedBlackGaussSeidel(Fields2D &fields, int nx, int ny, double coef,
                              int maxIters, double tol, double beta) {
  fields.Div();
  double res0 = 1.0;

  for (int it = 0; it < maxIters; ++it) {
    for (int color = 0; color < 2; ++color) {
OMP_PRAGMA(omp parallel for collapse(2))
for (int j = 0; j < ny; ++j) {
  for (int i = 0; i < nx; ++i) {
    if ((i + j) % 2 != color)
      continue;
    fields.p.Set(i, j, gsUpdate(fields, nx, ny, i, j, coef, beta));
  }
}
    }

    const double res = residualNorm(fields, nx, ny, coef, beta);
    //DBG_PRINTF("%f" ,res0!=0.0? res / res0:0);
    
    if (checkConvergence(res, res0, it, tol)) {
#ifndef NDEBUG
      std::cout << "  RedBlackGS converged in " << it + 1
                << " iters, rel.res = " << res / res0 << '\n';
#endif
      return;
    }
  }
#ifndef NDEBUG
  std::cout << "  RedBlackGS: reached maxIters = " << maxIters << '\n';
#endif
}

bool solveMICCG0(Fields2D &fields, double scale, int maxIters, double tol){
  fields.Div();

  const int nx = fields.nx;
  const int ny = fields.ny;
  const int N = nx * ny;

  // allocated once before and reused
  std::vector<double> Adiag(N, 0.0);

  //   Persistent : Adiag, precon, p, r, s
  //   Scratch    : z, q, tmp
  //
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (fields.Label(i, j) != Fields2D::FLUID)
        continue;
      const int id = nx * j + i;
      if (i - 1 >= 0 && fields.Label(i - 1, j) != Fields2D::SOLID)
        Adiag[id] += scale;
      if (i + 1 < nx && fields.Label(i + 1, j) != Fields2D::SOLID)
        Adiag[id] += scale;
      if (j - 1 >= 0 && fields.Label(i, j - 1) != Fields2D::SOLID)
        Adiag[id] += scale;
      if (j + 1 < ny && fields.Label(i, j + 1) != Fields2D::SOLID)
        Adiag[id] += scale;
    }
  }

  // residual r= b = -div at initial guess
  std::vector<double> r(N, 0.0);

  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      if (fields.Label(i, j) == Fields2D::FLUID)
        r[nx * j + i] = -static_cast<double>(fields.div.Get(i, j));

  double r0_inf = 0.0;

  OMP_PRAGMA(omp parallel for reduction(max:r0_inf))
  for (int k = 0; k < N; ++k)
    r0_inf = std::max(r0_inf, std::abs(r[k]));
  if (r0_inf < 1e-30) {
#ifndef NDEBUG
    std::cout << "  MICCG0: initial residual negligible, skipping solve.\n";
#endif
    return true;
  }

  std::vector<double> precon(N, 0.0);
  buildPrecon(fields, nx, ny, scale, Adiag, precon);

#ifndef NDEBUG
  {
    int zeroDiag = 0;
    for (int jj = 0; jj < ny; ++jj)
      for (int ii = 0; ii < nx; ++ii)
        if (fields.Label(ii, jj) == Fields2D::FLUID &&
            Adiag[nx * jj + ii] == 0.0)
          ++zeroDiag;
    if (zeroDiag > 0)
      std::cout << "  MICCG0 WARNING: " << zeroDiag
                << " fluid cells have Adiag=0\n";
  }
#endif

  // one allocation, reused every iteration
  std::vector<double> p(N, 0.0);
  std::vector<double> z(N, 0.0);
  std::vector<double> s(N, 0.0);
  std::vector<double> tmp(N, 0.0);
  std::vector<double> q(N, 0.0);

  // PCG initialization
  // CF fig 5.6
  applyPrecon(fields, nx, ny, scale, precon, r, z, q);
  s = z;

  double sigma_pcg = 0.0;
  OMP_PRAGMA(omp parallel for reduction(+:sigma_pcg))
  for (int k = 0; k < N; ++k)
    sigma_pcg += z[k] * r[k];

#ifndef NDEBUG
  if (sigma_pcg <= 0.0)
    std::cout << "  MICCG0 WARNING: sigma_pcg = " << sigma_pcg
              << " (preconditioner not SPD)\n";
#endif
  // PCG
  bool converged = false;
  // loop until done or max iter exceed
  for (int it = 0; it < maxIters; ++it) {
    // set auxiliary vector z= applyA(s)
    applyA(fields, nx, ny, scale, Adiag, s, tmp);
    // alpha = sigma / dot(z,s)
    double dot_s_tmp = 0.0;
    OMP_PRAGMA(omp parallel for reduction(+:dot_s_tmp))
    for (int k = 0; k < N; ++k)
      dot_s_tmp += tmp[k] * s[k];
    // round off
    if (std::abs(dot_s_tmp) < 1e-30)
      break;

    const double alpha = sigma_pcg / dot_s_tmp;

    double r_inf = 0.0;
    // fused for a single loop
  OMP_PRAGMA(omp parallel for reduction(max:r_inf))
  for (int k = 0; k < N; ++k) {
    p[k] += alpha * s[k];
    r[k] -= alpha * tmp[k];
    r_inf = std::max(r_inf, std::abs(r[k]));
  }

  // if max r <= tol then return p
  if (r_inf / r0_inf <= tol) {
#ifndef NDEBUG
    std::cout << "  MICCG0 converged in " << it + 1
              << " iters, rel.res = " << r_inf / r0_inf << '\n';
#endif
    converged = true;
    break;
  }
  // otherwise failsafe
  if (it == maxIters - 1) {
#ifndef NDEBUG
    std::cout << "  MICCG0: reached maxIters = " << maxIters
              << ", rel.res = " << r_inf / r0_inf << '\n';
#endif
    break;
  }

  applyPrecon(fields, nx, ny, scale, precon, r, z, q);

  double sigma_new = 0.0;
    OMP_PRAGMA(omp parallel for reduction(+:sigma_new))
    for (int k = 0; k < N; ++k)
      sigma_new += z[k] * r[k];

    const double beta = sigma_new / sigma_pcg;
    sigma_pcg = sigma_new;
    OMP_PRAGMA(omp parallel for)
    for (int k = 0; k < N; ++k)
      s[k] = z[k] + beta * s[k];
  }
  // assumed p is good now
  // TODO fields.p can be changed with p[] if p is a grid2D
  // write solution back
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      if (fields.Label(i, j) == Fields2D::FLUID)
        fields.p.Set(i, j, static_cast<varType>(p[nx * j + i]));

  return converged;
}
