#pragma once
// utils.hpp  —  stencil helpers shared by RBGS.cpp and CG.cpp (CPU only).
// Do NOT include from .cu files.

#include "../Fields.hpp"
#include "../Precision.hpp"
#include <cassert>
#include <cmath>
#include <limits>

/* small reminder about the geometry
 *           +----------+
 *           |          |
 *           |  P       |
 *           |   i,j+1  |
 *           |          |
 *           |          |
 * +---------+-V--------+----------+
 * |         |  i-1,j   |          |
 * | P       |  P       |  P       |
 * |  i-1,j  |   i,j    |   i+1,j  |
 * |         U          U          |
 * |         |i-1,j-1   |i,j-1     |
 * +---------+-V--------+----------+
 *           |  i-1,j-1 |
 *           |          |
 *           |  P       |
 *           |   i,j-1  |
 *           |          |
 *           +----------+
 */
[[nodiscard]] inline varType neighbourSum(const Fields2D &f,
                                          const int i, const int j,
                                          [[maybe_unused]] varType beta) noexcept {
    assert(i >= 1 && i < f.p.nx - 1);
    assert(j >= 1 && j < f.p.ny - 1);

    const varType pC  = f.p.Get(i, j);
    assert(std::isfinite(pC));

    const labeltype cur = f.Label(i, j);
    varType sumP = static_cast<varType>(0);

    // Left (i-1, j)
    { const labeltype nb = f.Label(i - 1, j);
      if      (IS_SOLID(nb) || IS_BC_U(nb)) sumP += pC;
      else if (!IS_AIR(nb))                 sumP += f.p.Get(i - 1, j); }

    // Right (i+1, j)
    { const labeltype nb = f.Label(i + 1, j);
      if      (IS_SOLID(nb) || IS_BC_U(cur)) sumP += pC;
      else if (!IS_AIR(nb))                  sumP += f.p.Get(i + 1, j); }

    // Bottom (i, j-1)
    { const labeltype nb = f.Label(i, j - 1);
      if      (IS_SOLID(nb) || IS_BC_V(nb)) sumP += pC;
      else if (!IS_AIR(nb))                 sumP += f.p.Get(i, j - 1); }

    // Top (i, j+1)
    { const labeltype nb = f.Label(i, j + 1);
      if      (IS_SOLID(nb) || IS_BC_V(cur)) sumP += pC;
      else if (!IS_AIR(nb))                  sumP += f.p.Get(i, j + 1); }

    return sumP;
}

[[nodiscard]] inline varType gsUpdate(const Fields2D &f,
                                     const int i, const int j,
                                     const varType coef,
                                     const varType beta) noexcept {
    return (-coef * f.div.Get(i - 1, j - 1)
            + neighbourSum(f, i, j, beta)) / 4.0;
}

inline bool checkConvergence(varType res, varType &res0,
                              int it, varType tol) noexcept {
    if (it == 0) {
        res0 = res;
        return res0 < 1e-30;
    }
    return (res0 < 1e-30) || (res / res0 < tol);
}

inline void debugSolverConverged(const char *solverName, const int iterations,
                                 const double relRes) {
#ifdef NDEBUG
  (void)solverName;
  (void)iterations;
  (void)relRes;
#else
  DBG_PRINTF("%s converged in %d iters, rel.res = %.6g", solverName, iterations,
             relRes);
#endif
}

inline void debugSolverMaxIters(
    const char *solverName, const int maxIters,
    const double relRes = std::numeric_limits<double>::quiet_NaN()) {
#ifdef NDEBUG
  (void)solverName;
  (void)maxIters;
  (void)relRes;
#else
  if (std::isfinite(relRes)) {
    DBG_PRINTF("%s: reached maxIters = %d, rel.res = %.6g", solverName,
               maxIters, relRes);
  } else {
    DBG_PRINTF("%s: reached maxIters = %d", solverName, maxIters);
  }
#endif
}
