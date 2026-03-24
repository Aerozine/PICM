#pragma once
#include "Fields.hpp"

/**
 * @file IterativeMethods.hpp
 * @brief All iterative method implemented ( jacobi , RB-GS , CG )
 *         optimized for poisson problem
 *
 * All functions operate directly on a @c Fields2D instance: they read the
 * cell labels and divergence field, and write the resulting pressure back.
 * They have no dependency on any particular solver class.
 *
 * ### Scaling convention (Jacobi / GS / RBGS)
 * These three solvers expect
 * @code
 *   coef = rho * dx^2 / dt
 * @endcode
 * which is the natural right-hand-side scaling for the GS update formula.
 *
 * ### Scaling convention (MICCG0)
 * MICCG0 expects the complementary value
 * @code
 *   scale = dt / (rho * dx^2)
 * @endcode
 * With this convention the system is A·p = b where b = −div (raw negative
 * divergence, units 1/s), and the solution p is in Pascals, consistent with
 * the velocity correction
 * @code
 *   u -= (dt / (rho * dx)) * (p[i] - p[i-1])
 * @endcode
 */

/**
 * @brief Jacobi pressure solver (fully parallel, slower convergence).
 *
 * @param fields    Simulation fields @c div read, @c p updated in-place.
 * @param nx        Number of pressure cells in x.
 * @param ny        Number of pressure cells in y.
 * @param coef      rho * dx^2 / dt.
 * @param maxIters  Maximum number of iterations.
 * @param tol       Relative residual convergence threshold.
 */
void solveJacobi(Fields2D &fields, int nx, int ny, double coef, int maxIters,
                 double tol);

/**
 * @brief Gauss-Seidel pressure solver (sequential, faster convergence).
 *
 * @param fields    Simulation fields @c div read, @c p updated in-place.
 * @param nx        Number of pressure cells in x.
 * @param ny        Number of pressure cells in y.
 * @param coef      rho * dx^2 / dt.
 * @param maxIters  Maximum number of iterations.
 * @param tol       Relative residual convergence threshold.
 */
void solveGaussSeidel(Fields2D &fields, int nx, int ny, double coef,
                      int maxIters, double tol);

/**
 * @brief Red-Black Gauss-Seidel pressure solver (parallel + fast convergence).
 *
 * @param fields    Simulation fields @c div read, @c p updated in-place.
 * @param nx        Number of pressure cells in x.
 * @param ny        Number of pressure cells in y.
 * @param coef      rho * dx^2 / dt.
 * @param maxIters  Maximum number of iterations.
 * @param tol       Relative residual convergence threshold.
 */
void solveRedBlackGaussSeidel(Fields2D &fields, int nx, int ny, double coef,
                              int maxIters, double tol);

/**
 * @brief Solve the pressure Poisson equation using MICCG(0).
 * Modified Incomplete Cholesky Conjugate Gradient, Level 0
 * Convergence is measured with the infinity norm of the residual
 * @param fields    Simulation fields @c div read, @c p overwritten.
 * @param scale     dt / (rho * dx^2).
 * @param maxIters  Maximum number of PCG iterations.
 * @param tol       Relative residual convergence threshold.
 * @return          @c true if converged within @p maxIters.
 */
bool solveMICCG0(Fields2D &fields, double scale, int maxIters, double tol);
