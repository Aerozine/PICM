#pragma once

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
[[nodiscard]] inline varType neighbourSum(const Fields2D &f, const int i,
                                          const int j) noexcept {
  assert(i >= 1 && i < f.p.nx - 1);
  assert(j >= 1 && j < f.p.ny - 1);

  const varType pC = f.p.Get(i, j);
  assert(std::isfinite(pC));

  varType sumP = static_cast<varType>(0);

  // Left (i-1, j)
  {
    const labeltype nb = f.Label(i - 1, j);
    if (IS_SOLID(nb) || IS_BC_U(nb))
      sumP += pC;
    else if (IS_AIR(nb) && f.phi != nullptr)
      sumP += f.interface_u->Get(i - 1, j - 1);
    else if (!IS_AIR(nb))
      sumP += f.p.Get(i - 1, j);
  }

  // Right (i+1, j) — BC_U flag lives on the neighbour, not the current cell
  {
    const labeltype nb = f.Label(i + 1, j);
    if (IS_SOLID(nb) || IS_BC_U(nb))
      sumP += pC;
    else if (IS_AIR(nb) && f.phi != nullptr)
      sumP += f.interface_u->Get(i, j - 1);
    else if (!IS_AIR(nb))
      sumP += f.p.Get(i + 1, j);
  }

  // Bottom (i, j-1)
  {
    const labeltype nb = f.Label(i, j - 1);
    if (IS_SOLID(nb) || IS_BC_V(nb))
      sumP += pC;
    else if (IS_AIR(nb) && f.phi != nullptr)
      sumP += f.interface_v->Get(i - 1, j - 1);
    else if (!IS_AIR(nb))
      sumP += f.p.Get(i, j - 1);
  }

  // Top (i, j+1) — BC_V flag lives on the neighbour, not the current cell
  {
    const labeltype nb = f.Label(i, j + 1);
    if (IS_SOLID(nb) || IS_BC_V(nb))
      sumP += pC;
    else if (IS_AIR(nb) && f.phi != nullptr)
      sumP += f.interface_v->Get(i - 1, j);
    else if (!IS_AIR(nb))
      sumP += f.p.Get(i, j + 1);
  }

  return sumP;
}

[[nodiscard]] inline varType gsUpdate(const Fields2D &f, const int i,
                                      const int j,
                                      const varType coef) noexcept {
  return (-coef * f.div.Get(i - 1, j - 1) + neighbourSum(f, i, j)) / 4.0;
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
