
#include "../Precision.hpp"
#include "IterativeMethods.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

void solveRedBlackGaussSeidel(Fields2D &fields, int nx, int ny, varType coef,
                              int maxIters, varType tol, varType /*beta*/) {
  fields.Div();

  constexpr varType PI = 3.14159265358979323846;
  const int N_min = std::max(1, std::min(nx, ny));
  const varType omega = std::min(1.95, 2.0 / (1.0 + std::sin(PI / N_min)));

  const int pnx = fields.p.nx;
  const int pny = fields.p.ny;

  const varType norm = computeRhsNorm(fields, coef);

  // Each colour sweep is data-race-free : safe for OpenMP collapse.
  for (int it = 0; it < maxIters; ++it) {
    varType sumSq = 0.0;
    for (int colour = 0; colour < 2; ++colour) {
      OMP_PRAGMA(omp parallel for reduction(+:sumSq) schedule(static))
      // loop simplified to improve cache locality n branch prediction
      for (int j = 1; j < pny - 1; ++j) {
        const int i_start = 1 + ((j + colour + 1) & 1);
        for (int i = i_start; i < pnx - 1; i += 2) {
          const labeltype lbl = fields.Label(i, j);

          if (IS_SOLID(lbl) || IS_BC_P(lbl) || IS_AIR(lbl))
            continue;

          const varType p_old = fields.p.Get(i, j);
          const varType p_gs = gsUpdate(fields, i, j, coef);

          const varType p_new = p_old + omega * (p_gs - p_old);
          assert(std::isfinite(p_new));
          fields.p.Set(i, j, p_new);

          // True residual r = b + nS - 4*p_old  (= 4*(p_gs - p_old))
          const varType r = varType(4) * (p_gs - p_old);
          sumSq += r * r;
        }
      }
    }

    const varType residue = std::sqrt(sumSq);
    const varType condition = residue / norm;

    if (condition < tol) {
      debugSolverConverged("RedBlackGS", it + 1, condition);
      return;
    }
  }

  debugSolverMaxIters("RedBlackGS", maxIters);
}