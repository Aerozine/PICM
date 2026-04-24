#pragma once
// utils.hpp  —  stencil helpers shared by RBGS.cpp and CG.cpp (CPU only).
// Do NOT include from .cu files.

#include "../Fields.hpp"
#include "../Precision.hpp"
#include <cassert>
#include <cmath>

/* Stencil geometry (P-space)
 *
 *           +----------+
 *           |  P i,j+1 |
 * +---------+----V-----+----------+
 * | P i-1,j |  P i,j   | P i+1,j |
 *           U          U
 * +---------+----V-----+----------+
 *           |  P i,j-1 |
 *           +----------+
 *
 *  BC_U label lives on the face between (i-1,j) and (i,j)  → Label(i, j)
 *  BC_V label lives on the face between (i,j-1) and (i,j)  → Label(i, j)
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
