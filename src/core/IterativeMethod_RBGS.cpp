#include "Fields.hpp"
#include "Grid2D.hpp"
#include "IterativeMethods.hpp"
#include "Precision.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

//@todo I know this should be illegal
// but this is the simpliest way to incude code and nicely inside a c file
#include "IterativeMethods_utils.cpp"

void solveRedBlackGaussSeidel(Fields2D &fields, int nx, int ny, double coef,
                              int maxIters, double tol, double beta) {
  fields.Div();
  double res0 = 1.0;

  const int N = std::min(nx, ny);
  // @todo this is really overkill apollo used the same number
  // of decimal to reach the moon
  constexpr double pi = 3.14159265358979;
  const double omega = std::min(1.95, 2.0 / (1.0 + std::sin(pi / N)));

  for (int it = 0; it < maxIters; ++it) {
    double sumSq = 0.0;
    int count = 0;
    // inner domain nx ny from fields
    // @todo rewrite color as
    // n =0 , n+=2 , n<nmax
    // n =1 , n+=2 , n<nmax
    // for cache optim and simplification
    // ALL INDICIES IN P referential
    /*
    if (fields.Label(i + 1, j + 1) & Fields2D::SOLID ||
      fields.Label(i + 1, j + 1) & Fields2D::BC_P ||
      fields.Label(i + 1, j + 1) & Fields2D::AIR ||
      fields.Label(i + 1, j + 1) & Fields2D::IC_P) {
      @todo trying IC_P is idiomatic because it does a BC job
      does BCP is really correct ?
      */
    // in P system
    for (int color = 0; color < 2; ++color) {
OMP_PRAGMA(omp parallel for collapse(2) reduction(+:sumSq) reduction(+:count))
// for the inner solution ;
for (int j = 1; j < fields.p.ny - 1; ++j) {
  for (int i = 1; i < fields.p.nx - 1; ++i) {
    if ((i + j) % 2 != color)
      continue;
    labeltype current_label = fields.Label(i, j);
    if (IS_SOLID(current_label) || IS_BC_P(current_label) ||
        IS_AIR(current_label))
      continue;
    const double p_old = fields.p.Get(i, j);
    const auto sumP =
        neighbourSum(fields, i, j, beta);
    // i j express in p ref , for div which is nx , ny -1,-1
    // @todo always /4 but not always 4 neighbour right ?
    const double p_neighbour =
        (-coef * fields.div.Get(i - 1, j - 1) + sumP) / 4.0;
    const double p_new = p_old + omega * (p_neighbour - p_old);
    assert(std::isfinite(p_new));
    fields.p.Set(i, j, p_new);

    const double r = p_new - p_old;
    sumSq += r * r;
    ++count;
  }
}
    }
    // assert(sumSq > 0.0);

    const double res = (count > 0) ? std::sqrt(sumSq / count) : 0.0;
    if (checkConvergence(res, res0, it, tol)) {
      // @todo use DEBUG PRINT MACRO please
#ifndef NDEBUG
      std::cout << "  RedBlackGS converged in " << it + 1
                << " iters, rel.res = " << res / res0 << '\n';
      assert(std::isfinite(res / res0));
#endif
      return;
    }
  }
  // #ifndef NDEBUG
  std::cout << "  RedBlackGS: reached maxIters = " << maxIters << '\n';
  // #endif
}
