#pragma once
#include "../../core/Fields.hpp"
#include "../../core/OutputWriter.hpp"
#include "../../core/Parameters.hpp"
#include <memory>

/**
 * @file SemiLagrangian.hpp
 * @brief Semi-Lagrangian incompressible Navier-Stokes solver on a MAC grid.
 */

/**
 * @brief 2-D incompressible Navier-Stokes solver using a semi-Lagrangian
 *        advection scheme and a pressure-projection method.
 *
 * ### Algorithm — one time step
 * 1. **Project** (+MakeIncompressible): solve the pressure Poisson equation
 *    and correct velocities so that \f$\nabla \cdot \mathbf{u} \approx 0 \f$.
 * 2. **Advect**: trace departure points backward in time (RK2) and
 *    interpolate the velocity field at those points.
 */
class SemiLagrangian {
public:
  /**
   * @brief Construct the solver, initialise fields, and open output writers.
   * @param params Simulation parameters (non-owning reference, must outlive
   *               this object).
   */
  explicit SemiLagrangian(const Parameters &params);

  ~SemiLagrangian();

  SemiLagrangian(const SemiLagrangian &) = delete;
  SemiLagrangian &operator=(const SemiLagrangian &) = delete;

  /// @brief Run the full simulation loop (nt steps) and write output.
  void Run();

  /// @brief Advance the simulation by one time step.
  void Step();

  Fields2D &GetFields() { return *fields; } ///< Access fields (mutable).
  const Fields2D &GetFields() const {
    return *fields;
  } ///< Access fields (const).

private:
  const Parameters &params;

  // Cached scalars from params to avoid pointer chasing in hot loops.
  int nx, ny;
  varType dx, dy, dt;
  varType density;

  Fields2D *fields; ///< @todo Replace with std::unique_ptr<Fields2D>.

  // Output writers — null if the corresponding write_* flag is false.
  std::unique_ptr<OutputWriter> uWriter;
  std::unique_ptr<OutputWriter> vWriter;
  std::unique_ptr<OutputWriter> pWriter;
  std::unique_ptr<OutputWriter> divWriter;
  std::unique_ptr<OutputWriter> normVelocityWriter;

  /// @brief Construct the OutputWriters requested in @c params.
  void InitializeOutputWriters();

  /**
   * @brief Write all enabled fields at the current step if it falls on a
   *        sampling interval.
   * @param step Current time-step index (0-based).
   */
  void WriteOutput(int step) const;

  // Advection

  /**
   * @brief Advect u and v using a semi-Lagrangian (RK2 backward-trace +
   *        bilinear interpolation) scheme.
   */
  void Advect() const;

  /**
   * @brief Trace the departure point of a u-face at grid position (i, j)
   *        backward in time using RK2.
   *
   * The u-face is located at physical position (i.dx, (j+0.5).dy).
   *
   * @param[in]  i  Face x-index.
   * @param[in]  j  Face y-index.
   * @param[out] x  Physical x-coordinate of the departure point.
   * @param[out] y  Physical y-coordinate of the departure point.
   */
  void traceParticleU(int i, int j, varType &x, varType &y) const;

  /**
   * @brief Trace the departure point of a v-face at grid position (i, j)
   *        backward in time using RK2.
   *
   * The v-face is located at physical position ((i+0.5).dx, j.dy).
   *
   * @param[in]  i  Face x-index.
   * @param[in]  j  Face y-index.
   * @param[out] x  Physical x-coordinate of the departure point.
   * @param[out] y  Physical y-coordinate of the departure point.
   */
  void traceParticleV(int i, int j, varType &x, varType &y) const;

  /**
   * @brief Bilinearly interpolate the u field at physical position (x, y).
   * @param x Physical x-coordinate (clamped to the domain).
   * @param y Physical y-coordinate (clamped to the domain).
   * @return  Interpolated u value.
   */
  [[nodiscard]] varType interpolateU(varType x, varType y) const;

  /**
   * @brief Bilinearly interpolate the v field at physical position (x, y).
   * @param x Physical x-coordinate (clamped to the domain).
   * @param y Physical y-coordinate (clamped to the domain).
   * @return  Interpolated v value.
   */
  [[nodiscard]] varType interpolateV(varType x, varType y) const;

  /**
   * @brief Return both velocity components at physical position (x, y).
   * @param[in]  x Physical x-coordinate.
   * @param[in]  y Physical y-coordinate.
   * @param[out] u Interpolated u value.
   * @param[out] v Interpolated v value.
   */
  void getVelocity(varType x, varType y, varType &u, varType &v) const;

  // Projection
  /**
   * @brief Enforce \f$ \nabla \cdot \mathbf{u} = 0 \f$: solve pressure, then
   * correct velocities.
   */
  void MakeIncompressible();

  /**
   * @brief Dispatch to the pressure solver selected in @c params.
   * @param maxIters Maximum number of solver iterations.
   * @param tol      Relative residual convergence threshold.
   */
  void solvePressure(int maxIters, double tol);

  /**
   * @brief Apply the pressure gradient to correct face velocities.
   *
   * Implements the explicit update:
   * \f [ u^{n+1} = u^* - \frac{\Delta t}{\rho\,\Delta x}\,(p_i - p_{i-1}) \f]
   * Faces adjacent to SOLID cells are set to @c usolid instead.
   */
  void updateVelocities();

  /**
   * @brief Compute the RMS residual of the discrete Poisson equation.
   *
   * The residual at each FLUID cell is:
   * \f$ r_{ij} = -\text{coef}\cdot\text{div}_{ij}
   *              + \sum_{\text{nb}} p_{\text{nb}}
   *              - N\,p_{ij} \f$
   *
   * @param coef  Scaling coefficient \f$\rho\,\Delta x^2 / \Delta t \f$.
   * @return RMS residual over all FLUID cells (0 if none).
   */
  [[nodiscard]] double computeResidualNorm(varType coef) const;

  /**
   * @brief Compute the Gauss-Seidel update for cell (i, j).
   *
   * \f$ p^{\text{new}}_{ij} =
   *     \frac{-\text{coef}\cdot\text{div}_{ij} + \sum_{\text{nb}}
   * p_{\text{nb}}}{N} \f$
   *
   * @param i    Cell x-index.
   * @param j    Cell y-index.
   * @param coef Scaling coefficient.
   * @return     New pressure value, or NAN if the cell is not FLUID.
   */
  [[nodiscard]] double getUpdate(int i, int j, varType coef) const;

  /// @brief Jacobi pressure solver (fully parallel, slower convergence).
  void SolveJacobi(int maxIters, double tol);

  /// @brief Gauss-Seidel pressure solver (sequential, faster convergence).
  void SolveGaussSeidel(int maxIters, double tol);

  /// @brief Red-Black Gauss-Seidel pressure solver (parallel + fast
  /// convergence).
  void SolveRedBlackGaussSeidel(int maxIters, double tol);
};
