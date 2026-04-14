#include "Fields.hpp"
#include "Grid2D.hpp"
#include "IterativeMethods.hpp"
#include "Precision.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
// should be statitced
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

// IN P SPACE PLEASE
[[nodiscard]] inline varType neighbourSum(const Fields2D &f, int nx, int ny,
                                          int i, int j, double beta) {
  // i j in terms of pressure cells
  varType sumP = 0.0;
  const varType pC = f.p.Get(i, j);
  assert(std::isfinite(pC));
  assert(i >= 0 && i < nx);
  assert(j >= 0 && j < ny);

  labeltype current_cell = f.Label(i, j);
  // Left neighbour
  labeltype left_cell = f.Label(i - 1, j);
  // @todo to be checked correctly (first condition with the value
  if (IS_SOLID(left_cell)) {
    sumP += pC;
  } else if (IS_BC_U(left_cell)) {
    // in ghost = pC- beta*(-u_BC)
    //  u_bc is j-1 and left i-1
    sumP += pC + beta * f.u.Get(i - 1, j - 1);
  } else if (!IS_AIR(left_cell))
    sumP += f.p.Get(i - 1, j);

  // Right neigbour
  labeltype right_cell = f.Label(i + 1, j);
  if (IS_SOLID(right_cell)) {
    sumP += pC;
  } else if (IS_BC_U(current_cell)) {
    // out ghost= pC - beta*u_bc
    //  curent cell is i,j-1
    sumP += pC - beta * f.u.Get(i, j - 1);
  } else if (!IS_AIR(right_cell))
    sumP += f.p.Get(i + 1, j);

  // Bottom neighbour
  labeltype bottom_cell = f.Label(i, j - 1);
  if (IS_SOLID(bottom_cell)) {
    sumP += pC;
  } else if (IS_BC_V(bottom_cell)) {
    // in -> +
    sumP += pC + beta * f.v.Get(i - 1, j - 1);
  } else if (!IS_AIR(bottom_cell))
    sumP += f.p.Get(i, j - 1);

  // Top neighbour
  labeltype top_cell = f.Label(i, j + 1);
  if (IS_SOLID(top_cell)) {
    sumP += pC;
  } else if (IS_BC_V(current_cell)) {
    sumP += pC - beta * f.v.Get(i - 1, j);
  } else if (!IS_AIR(top_cell))
    sumP += f.p.Get(i, j + 1);

  return sumP;
}
// Relative convergence: records res0 on first call (it==0).
inline bool checkConvergence(double res, double &res0, int it, double tol) {
  if (it == 0) {
    res0 = res;
    return res0 < 1e-30;
  }
  return (res0 < 1e-30) ? true : (res / res0) < tol;
}
// Gauss-Seidel update for a single FLUID cell.
//  p_new = ( -coef * div_{ij} + sum p_nb ) / N_nb
[[nodiscard]] inline double gsUpdate(const Fields2D &f, int nx, int ny, int i,
                                     int j, double coef, double beta) {
  const auto sumP = neighbourSum(f, nx, ny, i, j, beta);
  return (-coef * f.div.Get(i, j) + sumP) / 4.0;
}
