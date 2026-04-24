
#include "IterativeMethods.hpp"
#include "Utils.hpp"
#include "../Precision.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

//@todo to be fixed with the actual config
void solveJacobi(Fields2D &fields, int /*nx*/, int /*ny*/, varType coef,
                 int maxIters, varType tol, varType beta) {
    fields.Div();

    Grid2D pNew(fields.p.nx, fields.p.ny);
    varType res0 = 1.0;

    for (int it = 0; it < maxIters; ++it) {
        OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
        for (int j = 1; j < fields.p.ny - 1; ++j)
            for (int i = 1; i < fields.p.nx - 1; ++i) {
                const labeltype lbl = fields.Label(i, j);
                if (IS_SOLID(lbl) || IS_BC_P(lbl) || IS_AIR(lbl)) continue;
                pNew.Set(i, j, static_cast<varType>(gsUpdate(fields, i, j, coef, beta)));
            }

        varType sumSq = 0.0;
        int    count = 0;

        OMP_PRAGMA(omp parallel for collapse(2) reduction(+:sumSq) reduction(+:count) schedule(static))
        for (int j = 1; j < fields.p.ny - 1; ++j)
            for (int i = 1; i < fields.p.nx - 1; ++i) {
                if (!IS_FLUID(fields.Label(i, j))) continue;
                const varType r = static_cast<varType>(pNew.Get(i, j))
                               - static_cast<varType>(fields.p.Get(i, j));
                fields.p.Set(i, j, pNew.Get(i, j));
                sumSq += r * r;
                ++count;
            }

        const varType res = (count > 0) ? std::sqrt(sumSq / count) : 0.0;
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

void solveGaussSeidel(Fields2D &fields, int /*nx*/, int /*ny*/, varType coef,
                      int maxIters, varType tol, varType beta) {
    fields.Div();
    varType res0 = 1.0;

    for (int it = 0; it < maxIters; ++it) {
        varType sumSq = 0.0;
        int    count = 0;

        for (int j = 1; j < fields.p.ny - 1; ++j) {
            for (int i = 1; i < fields.p.nx - 1; ++i) {
                const labeltype lbl = fields.Label(i, j);
                if (IS_SOLID(lbl) || IS_BC_P(lbl) || IS_AIR(lbl)) continue;
                const varType p_old = static_cast<varType>(fields.p.Get(i, j));
                const varType p_new = gsUpdate(fields, i, j, coef, beta);
                fields.p.Set(i, j, static_cast<varType>(p_new));
                const varType r = p_new - p_old;
                sumSq += r * r;
                ++count;
            }
        }

        const varType res = (count > 0) ? std::sqrt(sumSq / count) : 0.0;
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

