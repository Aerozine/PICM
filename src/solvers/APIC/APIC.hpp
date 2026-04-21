#pragma once
#include "../PIC/PIC.hpp"

// todo: If we want reseeding (von-karman), need to to adapt seeding.cpp
//       to initialize correctly cux, cuy, cvx and cvy

class APIC : public PIC {
public:
  APIC(Parameters &params);

protected:
  void ProjectGridOnParticles() override;
  void ScatterToGridAPIC(varType xg, varType yg, varType xp, varType yp, varType baseVal, varType cX, varType cY, varType faceOffsetX, varType faceOffsetY, Grid2D &sum, Grid2D &weight, int imax, int jmax);
  void ProjectParticleOnMAC(int idx) override;
};