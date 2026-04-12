#pragma once
#include "../../core/Fields.hpp"
#include "../../core/IterativeMethods.hpp"
#include "../../core/OutputWriter.hpp"
#include "../../core/Parameters.hpp"
#include "Project.hpp"
#include <memory>
class SemiLagrangian {
public:
  /**
   * @brief Construct the solver, initialise fields, and open output writers.
   * @param params Simulation parameters (non-owning reference, must outlive
   *               this object).
   */

  SemiLagrangian(Parameters &params);

  ~SemiLagrangian();

  // SemiLagrangian(const SemiLagrangian &) = delete;
  // SemiLagrangian &operator=(const SemiLagrangian &) = delete;

  /// @brief Run the full simulation loop (nt steps) and write output.
  void Run();

  /// @brief Advance the simulation by one time step.
  void Step();

  // Fields2D &GetFields() { return *fields; } ///< Access fields (mutable).
  // const Fields2D &GetFields() const {
  // return *fields;
  //} ///< Access fields (const).

private:
  Parameters &params;

  // Cached scalars from params to avoid pointer chasing in hot loops.
  int nx, ny;
  varType dx, dy, dt;
  varType density;

  Fields2D *fields; ///< @todo Replace with std::unique_ptr<Fields2D>?

  // @todo maybe useful to do a solver class , and herit code from it ?
  // Output writers — null if the corresponding write_* flag is false.
  std::unique_ptr<OutputWriter> uWriter;
  std::unique_ptr<OutputWriter> vWriter;
  std::unique_ptr<OutputWriter> pWriter;
  std::unique_ptr<OutputWriter> divWriter;
  std::unique_ptr<OutputWriter> normVelocityWriter;
  std::unique_ptr<OutputWriter> smokeWriter;
  // @todo do a IF DEBUG THEN
  std::unique_ptr<OutputWriter> labelWriter;
  /// @brief Construct the OutputWriters requested in @c params.
  void InitializeOutputWriters();

  /**
   * @brief Write all enabled fields at the current step if it falls on a
   *        sampling interval.
   * @param step Current time-step index (0-based).
   */
  void WriteOutput(int step) const;

  /**
   * @brief Advect u and v using a semi-Lagrangian (RK2 backward-trace +
   *        bilinear interpolation) scheme.
   */
  void Advect() const;

  /**
   * @brief Advect smokeMap using a semi-Lagrangian (RK2 backward-trace +
   *        bilinear interpolation) scheme.
   */
  void AdvectSmoke() const;

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
   * @brief Bilinearly interpolate the smoke field at physical position (x, y).
   * @param x Physical x-coordinate (clamped to the domain).
   * @param y Physical y-coordinate (clamped to the domain).
   * @return  Interpolated smoke value.
   */
  [[nodiscard]] varType interpolateSmoke(varType x, varType y) const;

  /**
   * @brief Return both velocity components at physical position (x, y).
   * @param[in]  x Physical x-coordinate.
   * @param[in]  y Physical y-coordinate.
   * @param[out] u Interpolated u value.
   * @param[out] v Interpolated v value.
   */
  void getVelocity(varType x, varType y, varType &u, varType &v) const;
};
