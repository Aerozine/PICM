#pragma once
#include "../../core/Particles.hpp"
#include "../Solver.hpp"
#include <memory>

class PIC : public Solver {
public:
  explicit PIC(Parameters &params);
  ~PIC() override;

  void Run() override;
  void Step() override;

protected:
  Particles *particles = nullptr;

  varType GetW();
  [[nodiscard]] varType hat(varType r) const;

  virtual void ProjectParticleOnMAC(int idx);
  void ProjectParticlesOnGrid();
  void ProjectBCOnParticles();
  virtual void ProjectGridOnParticles();

  void Advect() override; // particle advection — overrides Solver::Advect()

  [[nodiscard]] varType interpolateU(const Grid2D &g, varType x,
                                     varType y) const;
  [[nodiscard]] varType interpolateV(const Grid2D &g, varType x,
                                     varType y) const;

  void RefillParticles();
  void CountAliveParticles();
  varType rand01();

  void ApplyGravity() const;
  void UpdateCellState() const;

  void ScatterToGrid(varType xg, varType yg, varType val, Grid2D &sum,
                     Grid2D &weight, int imax, int jmax);

  void WriteOutput(int step) const;

private:
  std::unique_ptr<OutputWriter> particlesWriter;
};
