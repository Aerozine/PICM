#include "SemiLagrangian.hpp"
#include <cmath>
#include <iostream>

// Cell update
double SemiLagrangian::getUpdate(const int i, const int j,
                                 const varType coef) const {
  if (fields->Label(i, j) != Fields2D::FLUID)
    return NAN;

  // Accumulate neighbour pressures and count valid neighbours.
  double sumP = 0.0;
  int nb = 0;

  if (i + 1 < nx) { sumP += fields->p.Get(i + 1, j); ++nb; }
  if (i - 1 >= 0) { sumP += fields->p.Get(i - 1, j); ++nb; }
  if (j + 1 < ny) { sumP += fields->p.Get(i, j + 1); ++nb; }
  if (j - 1 >= 0) { sumP += fields->p.Get(i, j - 1); ++nb; }

  if (nb == 0)
    return NAN;

  // Gauss-Seidel update:
  //   p_new = ( -coef * div_{ij} + Σ p_nb ) / N
  return (-coef * fields->div.Get(i, j) + sumP) / nb;
}

// Residual norm

double SemiLagrangian::computeResidualNorm(const varType coef) const {
  // RMS of the discrete Poisson residual over all FLUID cells:
  //   r_{ij} = rhs_{ij} - (A·p)_{ij}
  //          = -coef·div_{ij}  -  (nb·p_{ij} - Σ p_nb)
  double sumSq = 0.0;
  int count = 0;

#pragma omp parallel for collapse(2) reduction(+ : sumSq) reduction(+ : count)
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (fields->Label(i, j) != Fields2D::FLUID)
        continue;

      double sumP = 0.0;
      int nb = 0;

      if (i + 1 < nx) { sumP += fields->p.Get(i + 1, j); ++nb; }
      if (i - 1 >= 0) { sumP += fields->p.Get(i - 1, j); ++nb; }
      if (j + 1 < ny) { sumP += fields->p.Get(i, j + 1); ++nb; }
      if (j - 1 >= 0) { sumP += fields->p.Get(i, j - 1); ++nb; }

      const double r =
          (-coef * fields->div.Get(i, j)) - (nb * fields->p.Get(i, j) - sumP);
      sumSq += r * r;
      ++count;
    }
  }

  return (count > 0) ? std::sqrt(sumSq / count) : 0.0;
}

// Convergence check

// Returns true when the solver should stop.
// On the first call (it == 0), records res0 as the reference residual so that
// all subsequent checks use a *relative* criterion: ||r_k|| / ||r_0|| < tol.
static bool checkConvergence(const double res, double &res0, const int it,
                             const double tol) {
  if (it == 0) {
    res0 = res;
    return (res0 < 1e-30); // already converged if initial residual is tiny
  }
  return (res / res0) < tol;
}

// Jacobi

void SemiLagrangian::SolveJacobi(int maxIters, double tol) {
  const varType coef = density * dx * dx / dt;
  fields->Div();

  // Jacobi requires a separate buffer because all reads must use the
  // previous-iteration values.
  Grid2D pNew(nx, ny);
  double res0 = 1.0;

  for (int it = 0; it < maxIters; ++it) {

#pragma omp parallel for collapse(2)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        pNew.Set(i, j, getUpdate(i, j, coef));

#pragma omp parallel for collapse(2)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        if (fields->Label(i, j) == Fields2D::FLUID)
          fields->p.Set(i, j, pNew.Get(i, j));

    const double res = computeResidualNorm(coef);
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

// Gauss-Seidel

void SemiLagrangian::SolveGaussSeidel(int maxIters, double tol) {
  const varType coef = density * dx * dx / dt;
  fields->Div();

  double res0 = 1.0;

  for (int it = 0; it < maxIters; ++it) {
    // Sequential sweep — each cell sees the latest neighbour values.
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const double newVal = getUpdate(i, j, coef);
        if (!std::isnan(newVal))
          fields->p.Set(i, j, newVal);
      }

    const double res = computeResidualNorm(coef);
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

// Red-Black Gauss-Seidel

void SemiLagrangian::SolveRedBlackGaussSeidel(int maxIters, double tol) {
  const varType coef = density * dx * dx / dt;
  fields->Div();

  double res0 = 1.0;

  for (int it = 0; it < maxIters; ++it) {
    // Two-colour decomposition: "red" cells (i+j even) and "black" cells
    // (i+j odd). Each colour's cells are independent of one another, so
    // the inner loop can be parallelised without data races.
    for (int color = 0; color < 2; ++color) {
#pragma omp parallel for collapse(2)
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          if ((i + j) % 2 != color)
            continue;
          const double newVal = getUpdate(i, j, coef);
          if (!std::isnan(newVal))
            fields->p.Set(i, j, newVal);
        }
      }
    }

    const double res = computeResidualNorm(coef);
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
