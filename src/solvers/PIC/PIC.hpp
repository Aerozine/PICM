#pragma once
#include "../../core/Fields.hpp"
#include "../../core/OutputWriter.hpp"
#include "../../core/Parameters.hpp"
#include "../../core/Particles.hpp"
#include <memory>

/**
 * @file PIC.hpp
 * @brief Particle-In-Cell solver for 2-D incompressible flow.
 *
 * ### Algorithm one time step
 * 1. P2G  transfer particle velocities to the MAC grid (hat kernel).
 * 2. Project solve pressure Poisson, correct face velocities
 * 3. G2P  interpolate corrected grid velocities back to particles.
 * 4. Advect move particles with RK2; kill those leaving the domain or
 *             entering a solid, register their slots in the free-list.
 * 5. Reseed repopulate underpopulated fluid cells from the free-list.
 *
 * ### Particle budget
 * The particle array is allocated with CAPACITY_FACTOR × nx×ny×ppcx×ppcy
 * slots.  The extra slots start dead and are pushed onto the free-list so
 * the reseeder can inject particles at an open inlet without waiting for
 * particles to die at the outlet.
 */
class PIC {
public:
  explicit PIC(const Parameters &params);
  ~PIC();

  PIC(const PIC &) = delete;
  PIC &operator=(const PIC &) = delete;

  void Run();
  void Step();

  Fields2D &GetFields() { return *fields; }
  const Fields2D &GetFields() const { return *fields; }

private:
  const Parameters &params;

  int nx, ny;
  varType dx, dy, dt;
  varType density;

  Fields2D *fields;
  Particles *particles;
  ParticleSlots *deadSlots;

  std::unique_ptr<OutputWriter> uWriter;
  std::unique_ptr<OutputWriter> vWriter;
  std::unique_ptr<OutputWriter> pWriter;
  std::unique_ptr<OutputWriter> divWriter;
  std::unique_ptr<OutputWriter> normVelocityWriter;
  std::unique_ptr<OutputWriter> smokeWriter;
  std::unique_ptr<OutputWriter> particlesWriter;

  void InitializeOutputWriters();
  void WriteOutput(int step) const;

  /**
   * @brief Build the initial free-list of dead particle slots.
   *
   * Kills grid-slot particles whose home cell is solid, and pushes all
   * extra capacity slots (indices >= nx*ny*ppcx*ppcy) onto the free-list.
   * Must be called after applyToFields() and InitParticleGrid().
   */
  void InitFreeSlots();

  // P2G
  varType GetW();
  varType hat(varType r);
  void ProjectOneParticleOnMAC(varType x, varType y, varType up, varType vp);
  void ProjectParticlesOnGrid(std::string kernel);

  // G2P
  void ProjectGridOnParticles();

  // Advection
  void AdvectParticles();
  [[nodiscard]] varType interpolateU(varType x, varType y) const;
  [[nodiscard]] varType interpolateV(varType x, varType y) const;
  void getVelocity(varType x, varType y, varType &u, varType &v) const;

  // Reseeding
  void RefillParticles();
  void CountAliveParticles();
  varType rand01();

  // Pressure projection
  void MakeIncompressible();
  void solvePressure(int maxIters, double tol);
  void updateVelocities();
  [[nodiscard]] double computeResidualNorm(varType coef) const;
  [[nodiscard]] double getUpdate(int i, int j, varType coef) const;
  void SolveJacobi(int maxIters, double tol);
  void SolveGaussSeidel(int maxIters, double tol);
  void SolveRedBlackGaussSeidel(int maxIters, double tol);
};
