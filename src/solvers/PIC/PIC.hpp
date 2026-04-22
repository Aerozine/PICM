#pragma once
#include "../../core/Particles.hpp"
#include "../Solver.hpp"
#include <memory>

class PIC : public Solver {
public:
  explicit PIC(Parameters &params);

  void Run() override;
  void Step() override;

protected:
  //Particles *particles = nullptr;
  std::unique_ptr<Particles> particles;
  varType GetW();

  varType hat(varType r) const;
  varType dhat(varType r) const;

  varType hatTrilinear(varType r) const;
  varType hatQuadraticBSpline(varType r) const;

  varType dhatTrilinear(varType r) const;
  varType dhatQuadraticBSpline(varType r) const;

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

  void ApplyGravity() const;
  void UpdateCellState() const;

  virtual void ScatterToGrid(varType xg, varType yg, varType val, Grid2D &sum,
                     Grid2D &weight, int imax, int jmax);

  void WriteOutput(int step) const;

private:
  std::unique_ptr<OutputWriter> particlesWriter;
  std::unique_ptr<OutputWriter> countAliveParticles_writer ;
};
