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
inline varType GetW() {
  return static_cast<varType>(particles->ppcx * particles->ppcy);
}
inline varType hat(varType r) const {
#ifdef USE_SPEED
    if (r >= varType(0) && r <= varType(1))
        return varType(1) - r;
    if (r >= varType(-1) && r < varType(0))
        return varType(1) + r;
    return varType(0);
#else
    if (varType(-1.5) <= r && r < varType(-0.5))
        return varType(0.5) * (r + varType(3.0 / 2.0)) * (r + varType(1.5));
    if (varType(-0.5) <= r && r < varType(0.5))
        return varType(0.75) - r * r;
    if (varType(0.5) <= r && r < varType(1.5))
        return varType(0.5) * (varType(1.5) - r) * (varType(1.5) - r);
    return varType(0);
#endif
}

inline varType dhat(varType r) const {
#ifdef USE_SPEED
    if (r > varType(0) && r < varType(1))
        return -varType(1);
    if (r > varType(-1) && r < varType(0))
      return varType(1);
    return varType(0);
#else
    if (varType(-1.5) <= r && r < varType(-0.5))
        return r + static_cast<varType>(1.5);
    if (varType(-0.5) <= r && r < varType(0.5))
        return -static_cast<varType>(2.0) * r;
    if (varType(0.5) <= r && r < varType(1.5))
        return r - static_cast<varType>(1.5);
    return static_cast<varType>(0);
#endif
}


  virtual void ProjectParticleOnMAC(int idx);
  void ProjectParticlesOnGrid();
  void ProjectBCOnParticles();
  virtual void ProjectGridOnParticles();

  void Advect() override; // particle advection — overrides Solver::Advect()

  void RefillParticles();

  void UpdateCellState() const;

  virtual void ScatterToGrid(varType xg, varType yg, varType val, Grid2D &sum,
                     Grid2D &weight, int imax, int jmax);

  void WriteOutput(int step) const;

private:
  std::unique_ptr<OutputWriter> particlesWriter;
  std::unique_ptr<OutputWriter> countAliveParticles_writer ;
};
