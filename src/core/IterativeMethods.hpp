#pragma once
#include "Fields.hpp"

/**
 * @brief Jacobi pressure solver (fully parallel, slower convergence).
 *
 * @param fields    Simulation fields @c div read, @c p updated in-place.
 * @param nx        Number of pressure cells in x.
 * @param ny        Number of pressure cells in y.
 * @param coef      rho * dx^2 / dt.
 * @param maxIters  Maximum number of iterations.
 * @param tol       Relative residual convergence threshold.
 * @param beta
 */
void solveJacobi(Fields2D &fields, int nx, int ny, double coef, int maxIters,
                 double tol, double beta);

/**
 * @brief Gauss-Seidel pressure solver (sequential, faster convergence).
 *
 * @param fields    Simulation fields @c div read, @c p updated in-place.
 * @param nx        Number of pressure cells in x.
 * @param ny        Number of pressure cells in y.
 * @param coef      rho * dx^2 / dt.
 * @param maxIters  Maximum number of iterations.
 * @param tol       Relative residual convergence threshold.
 * @param beta      Beta parameter for solid boundary conditions.
 */
void solveGaussSeidel(Fields2D &fields, int nx, int ny, double coef,
                      int maxIters, double tol, double beta);

/**
 * @brief Red-Black Gauss-Seidel pressure solver (parallel + fast convergence).
 *
 * @param fields    Simulation fields @c div read, @c p updated in-place.
 * @param nx        Number of pressure cells in x.
 * @param ny        Number of pressure cells in y.
 * @param coef      rho * dx^2 / dt.
 * @param maxIters  Maximum number of iterations.
 * @param tol       Relative residual convergence threshold.
 * @param beta      Beta parameter for solid boundary conditions.
 */
void solveRedBlackGaussSeidel(Fields2D &fields, int nx, int ny, double coef,
                              int maxIters, double tol, double beta);

/**
 * @brief Solve the pressure Poisson equation using MICCG(0).
 * @param fields    Simulation fields @c div read, @c p overwritten.
 * @param scale     dt / (rho * dx^2).
 * @param maxIters  Maximum number of PCG iterations.
 * @param tol       Relative residual convergence threshold.
 * @return          @c true if converged within @p maxIters.
 */
bool solveMICCG0(Fields2D &fields, double scale, int maxIters, double tol);
/**
 * @brief Unpreconditioned Conjugate Gradient pressure solver.
 *
 * Derives A from the same RBGS stencil:
 *   (Ap)(i,j) = 4*p(i,j) - neighbourSum(i,j)
 *   b(i,j)    = -coef * div(i-1, j-1)
 *
 * @param fields   div read, p overwritten (fluid cells only).
 * @param coef     rho * dx² / dt   (same as RBGS).
 * @param beta     rho * dx  / dt   (kept for API consistency, unused).
 * @param maxIters Maximum CG iterations.
 * @param tol      Convergence: ||r||_inf / ||b||_inf < tol.
 * @return         true if converged within maxIters.
 */
// bool solveCG(Fields2D &fields, double coef, double beta,
// int maxIters, double tol);
