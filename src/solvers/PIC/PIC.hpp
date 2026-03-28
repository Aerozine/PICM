#pragma once
#include "../../core/Fields.hpp"
#include "../../core/OutputWriter.hpp"
#include "../../core/Parameters.hpp"
#include "../../core/Particles.hpp"
#include "../../solvers/SemiLagrangian/Project.hpp"
#include <memory>

/**
 * @file PIC.hpp
 * @brief Particle-In-Cell solver for 2-D incompressible flow.
 */
class PIC {
public:

  PIC(const Parameters &params);

  ~PIC();

  //PIC(const PIC &) = delete;
  //PIC &operator=(const PIC &) = delete;

  void Run();
  void Step();

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
};
