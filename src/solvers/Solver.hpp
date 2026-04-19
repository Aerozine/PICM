#pragma once
#include "../core/Fields.hpp"
#include "../core/IterativeMethods.hpp"
#include "../core/OutputWriter.hpp"
#include "../core/Parameters.hpp"
#include <memory>

class Solver {
public:
  explicit Solver(Parameters &params);
  virtual ~Solver();
  // basicaly a if anyone try to copy or assign , bang error
  Solver(const Solver &) = delete;
  Solver &operator=(const Solver &) = delete;

  virtual void Run() = 0;
  virtual void Step() = 0;

protected:
  Parameters &params;

  // Cached scalars from params to avoid pointer chasing in hot loops.
  int nx, ny;
  varType dx, dy, dt;
  varType density;

  Fields2D *fields = nullptr; ///< @todo Replace with std::unique_ptr<Fields2D>?

  // @todo maybe useful to do a solver class , and herit code from it ?
  // Output writers — null if the corresponding write_* flag is false.
  std::unique_ptr<OutputWriter> uWriter;
  std::unique_ptr<OutputWriter> vWriter;
  std::unique_ptr<OutputWriter> pWriter;
  std::unique_ptr<OutputWriter> divWriter;
  std::unique_ptr<OutputWriter> normVelocityWriter;
  std::unique_ptr<OutputWriter> countAliveParticles;
  std::unique_ptr<OutputWriter> smokeWriter;
  // @todo do a IF DEBUG THEN
  std::unique_ptr<OutputWriter> labelWriter;
  // Shared time-step loop used by Run() in every subclass.
  void RunLoop(int reportEvery);
  /// @brief Construct the OutputWriters requested in @c params.
  void InitializeOutputWriters();
  // @todo maybe put it in-another file or update vel + solve pressure in
  // another scope
  void MakeIncompressible(const Parameters &p, Fields2D &f);

  void updateVelocities(const Parameters &p, Fields2D &f);

  void solvePressure(const Parameters &p, Fields2D &f);
  /**
   * @brief Write all enabled fields at the current step if it falls on a
   *        sampling interval.
   * @param step Current time-step index (0-based).
   */
  virtual void WriteOutput(int step) const;

  /**
   * @brief Advect u and v using a semi-Lagrangian (RK2 backward-trace +
   *        bilinear interpolation) scheme.
   */
  virtual void Advect();

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
