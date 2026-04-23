#pragma once
#include "../core/Fields.hpp"
#include "../core/IterativeMethod/IterativeMethods.hpp"
#include "../core/OutputWriter.hpp"
#include "../core/Parameters.hpp"
#include <memory>

class Solver {
public:
  explicit Solver(Parameters &params);
  //virtual ~Solver();
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

  //Fields2D *fields = nullptr; ///< @todo Replace with std::unique_ptr<Fields2D>?
  std::unique_ptr<Fields2D> fields;
  // Output writers — null if the corresponding write_* flag is false.
  std::unique_ptr<OutputWriter> uWriter;
  std::unique_ptr<OutputWriter> vWriter;
  std::unique_ptr<OutputWriter> pWriter;
  std::unique_ptr<OutputWriter> divWriter;
  std::unique_ptr<OutputWriter> normVelocityWriter;
  std::unique_ptr<OutputWriter> countAliveParticles;
  std::unique_ptr<OutputWriter> labelWriter;
  // Shared time-step loop used by Run() in every subclass.
  void RunLoop(int reportEvery);

  /// @brief Construct the OutputWriters requested in @c params.
  void InitializeOutputWriters();

  void MakeIncompressible(const Parameters &p, Fields2D &f);

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


void traceParticle(const int i, const int j, varType &x, varType &y,
                           uint8_t field) const;
  /**
   * @brief Return both velocity components at physical position (x, y).
   * @param[in]  x Physical x-coordinate.
   * @param[in]  y Physical y-coordinate.
   * @param[out] u Interpolated u value.
   * @param[out] v Interpolated v value.
   */
  inline void getVelocity(varType x, varType y, varType &u, varType &v) const {
    u = fields->u.interpolate<0>(x, y, dx, dy);
    v = fields->v.interpolate<1>(x, y, dx, dy);
  };
};
