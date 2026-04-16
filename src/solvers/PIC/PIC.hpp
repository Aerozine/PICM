#pragma once
#include "../../core/Fields.hpp"
#include "../../core/OutputWriter.hpp"
#include "../../core/Parameters.hpp"
#include "../../core/Particles.hpp"
#include "../../solvers/SemiLagrangian/Project.hpp"
#include <memory>

class PIC {
public:
  PIC(const Parameters &params);
  virtual ~PIC();
  void Run();
  virtual void Step();

protected:
  const Parameters &params;

  int nx, ny;
  varType dx, dy, dt;
  varType density;
  Fields2D *fields;
  Particles *particles;

  // particles to grid projection
  varType GetW();
  varType hat(varType r) const;
  virtual void ProjectParticleOnMAC(int idx);
  void ProjectParticlesOnGrid(std::string kernel);
  void ProjectBCOnParticles();

  // grid to particles projection
  virtual void ProjectGridOnParticles();

  // advection
  void AdvectParticles();
  [[nodiscard]] varType interpolateU(const Grid2D& g, varType x, varType y) const;
  [[nodiscard]] varType interpolateV(const Grid2D& g, varType x, varType y) const;

  // reseeding
  void RefillParticles();
  void CountAliveParticles();
  varType rand01();

  // gravity
  void ApplyGravity();

  // cellstate
  void UpdateCellState();

  void ScatterToGrid(varType xg, varType yg, varType val, Grid2D &sum,
                        Grid2D &weight, int imax, int jmax);

private:
  std::unique_ptr<OutputWriter> uWriter;
  std::unique_ptr<OutputWriter> vWriter;
  std::unique_ptr<OutputWriter> pWriter;
  std::unique_ptr<OutputWriter> divWriter;
  std::unique_ptr<OutputWriter> normVelocityWriter;
  std::unique_ptr<OutputWriter> smokeWriter;
  std::unique_ptr<OutputWriter> particlesWriter;

  std::unique_ptr<OutputWriter> labelWriter;

  void InitializeOutputWriters();
  void WriteOutput(int step) const;
};