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

  // PIC(const PIC &) = delete;
  // PIC &operator=(const PIC &) = delete;

  void Run();
  void Step();

private:
  const Parameters &params;

  int nx, ny;
  varType dx, dy, dt;
  varType density;
  Fields2D *fields;
  Particles *particles;
  // @todo maybe useful to do a solver class , and herit code from it ?
  std::unique_ptr<OutputWriter> uWriter;
  std::unique_ptr<OutputWriter> vWriter;
  std::unique_ptr<OutputWriter> pWriter;
  std::unique_ptr<OutputWriter> divWriter;
  std::unique_ptr<OutputWriter> normVelocityWriter;
  std::unique_ptr<OutputWriter> smokeWriter;
  std::unique_ptr<OutputWriter> particlesWriter;

  // @todo do a IF DEBUG THEN
  std::unique_ptr<OutputWriter> labelWriter;

  void InitializeOutputWriters();
  void WriteOutput(int step) const;

  // P2G
  varType GetW();
  varType hat(varType r);
  void ProjectParticleOnMAC(varType x, varType y, varType up, varType vp);
  void ProjectParticlesOnGrid(std::string kernel);

  // G2P
  void ProjectGridOnParticles();

  // Advection
  void AdvectParticles();
  [[nodiscard]] varType interpolateU(varType x, varType y) const;
  [[nodiscard]] varType interpolateV(varType x, varType y) const;

  // Reseeding
  void RefillParticles();
  void CountAliveParticles();
  varType rand01();

  // Gravity
  void ApplyGravity();

  // CellState
  void UpdateCellState();

  void ScatterToGrid(varType xg, varType yg, varType val, Grid2D &sum,
                     Grid2D &weight, int imax, int jmax);
};