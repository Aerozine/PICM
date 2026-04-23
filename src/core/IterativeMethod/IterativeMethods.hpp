#pragma once
#include "../Fields.hpp"
/**
 * @brief Jacobi pressure solver (fully parallel, slower convergence).
 * @param fields   div read, p updated in-place.
 * @param nx       Number of interior pressure cells in x.
 * @param ny       Number of interior pressure cells in y.
 * @param coef     rho * dx² / dt.
 * @param maxIters Maximum iterations.
 * @param tol      Relative residual convergence threshold.
 * @param beta     rho * dx / dt  (solid BC correction term).
 */
void solveJacobi(Fields2D &fields, int nx, int ny, varType coef,
                 int maxIters, varType tol, varType beta);

/**
 * @brief Gauss-Seidel pressure solver (sequential, faster convergence).
 */
void solveGaussSeidel(Fields2D &fields, int nx, int ny, varType coef,
                      int maxIters, varType tol, varType beta);

/**
 * @brief Red-Black Gauss-Seidel pressure solver (OpenMP parallel).
 */
void solveRedBlackGaussSeidel(Fields2D &fields, int nx, int ny, varType coef,
                              int maxIters, varType tol, varType beta);

/**
 * @brief Modified Incomplete Cholesky CG (MICCG0) pressure solver.
 * @param scale  dt / (rho * dx²).
 * @return true if converged within maxIters.
 */

bool solveMICCG0(Fields2D &fields, varType scale, int maxIters, varType tol);

/**
 * @brief Unpreconditioned Conjugate Gradient pressure solver (CPU).
 * @param coef     rho * dx² / dt.
 * @param beta     Kept for API consistency, unused by CG.
 * @return true if converged within maxIters.
 */
bool solveCG(Fields2D &fields, varType coef, varType beta,
             int maxIters, varType tol);

#ifdef USE_CUDA

/**
 * @brief Red-Black Gauss-Seidel pressure solver on the GPU (CUDA).
 *
 * Drops in as a replacement for solveRedBlackGaussSeidel().
 * The function signature is deliberately identical to the CPU version.
 */
bool solveRedBlackGaussSeidel_GPU(Fields2D &fields, int nx, int ny,
                                  varType coef, int maxIters, varType tol,
                                  varType beta);

/**
 * @brief Conjugate Gradient pressure solver on the GPU (CUDA).
 *
 * Drops in as a replacement for solveCG().
 */
bool solveCG_GPU(Fields2D &fields, varType coef, varType beta,
                 int maxIters, varType tol);

#endif // USE_CUDA
