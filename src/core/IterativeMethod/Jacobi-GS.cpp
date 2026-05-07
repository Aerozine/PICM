
#include "../Precision.hpp"
#include "IterativeMethods.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

varType computeRhsNorm(Fields2D &fields, varType coef) {
  const int pnx = fields.p.nx;
  const int pny = fields.p.ny;

  varType norm = 0.0;

  OMP_PRAGMA(omp parallel for collapse(2) reduction(+:norm) schedule(static))
  for (int j = 1; j < pny - 1; ++j) {
    for (int i = 1; i < pnx - 1; ++i) {
      const labeltype lbl = fields.Label(i, j);
      if (IS_SOLID(lbl) || IS_BC_P(lbl) || IS_AIR(lbl))
        continue;

      const varType b = coef * fields.div.Get(i - 1, j - 1);
      norm += b * b;
    }
  }

  return std::sqrt(norm) + REAL_EPSILON;
}

void solveJacobi(Fields2D &fields, int /*nx*/, int /*ny*/, varType coef,
                 int maxIters, varType tol, varType /*beta*/) {
  fields.Div();

  const int pnx = fields.p.nx;
  const int pny = fields.p.ny;

  Grid2D pNew(pnx, pny);
  const varType norm = computeRhsNorm(fields, coef);

  for (int it = 0; it < maxIters; ++it) {
    varType sumSq = 0.0;

    OMP_PRAGMA(omp parallel for collapse(2) reduction(+:sumSq) schedule(static))
    for (int j = 1; j < pny - 1; ++j) {
      for (int i = 1; i < pnx - 1; ++i) {
        const labeltype lbl = fields.Label(i, j);
        if (IS_SOLID(lbl) || IS_BC_P(lbl) || IS_AIR(lbl))
          continue;

        const varType p_old = fields.p.Get(i, j);
        const varType p_gs = gsUpdate(fields, i, j, coef);

        const varType r = varType(4) * (p_gs - p_old);
        sumSq += r * r;

        pNew.Set(i, j, static_cast<varType>(p_gs));
      }
    }

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 1; j < pny - 1; ++j) {
      for (int i = 1; i < pnx - 1; ++i) {
        const labeltype lbl = fields.Label(i, j);
        if (IS_SOLID(lbl) || IS_BC_P(lbl) || IS_AIR(lbl))
          continue;

        fields.p.Set(i, j, pNew.Get(i, j));
      }
    }

    const varType residue = std::sqrt(sumSq);
    const varType condition = residue / norm;

    if (condition < tol) {
      debugSolverConverged("Jacobi", it + 1, condition);
      return;
    }
  }

  debugSolverMaxIters("Jacobi", maxIters);
}

void solveGaussSeidel(Fields2D &fields, int /*nx*/, int /*ny*/, varType coef,
                      int maxIters, varType tol, varType /*beta*/) {
  fields.Div();

  const int pnx = fields.p.nx;
  const int pny = fields.p.ny;
  const varType norm = computeRhsNorm(fields, coef);

  for (int it = 0; it < maxIters; ++it) {
    varType sumSq = 0.0;

    for (int j = 1; j < pny - 1; ++j) {
      for (int i = 1; i < pnx - 1; ++i) {
        const labeltype lbl = fields.Label(i, j);
        if (IS_SOLID(lbl) || IS_BC_P(lbl) || IS_AIR(lbl))
          continue;

        const varType p_old = fields.p.Get(i, j);
        const varType p_gs = gsUpdate(fields, i, j, coef);

        const varType r = varType(4) * (p_gs - p_old);
        sumSq += r * r;

        assert(std::isfinite(p_gs));
        fields.p.Set(i, j, static_cast<varType>(p_gs));
      }
    }

    const varType residue = std::sqrt(sumSq);
    const varType condition = residue / norm;

    if (condition < tol) {
      debugSolverConverged("GaussSeidel", it + 1, condition);
      return;
    }
  }

  debugSolverMaxIters("GaussSeidel", maxIters);
}
