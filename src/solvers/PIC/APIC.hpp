#pragma once
#include "PIC.hpp"

// todo: If we want reseeding (von-karman), need to to adapt seeding.cpp
//       to initialize correctly cux, cuy, cvx and cvy

class APIC: public PIC {
public:
  APIC(const Parameters &params);

  void ProjectGridOnParticles() override;
  Vec2 gradWeightU(int i, int j, varType xp, varType yp) const;
  Vec2 gradWeightV(int i, int j, varType xp, varType yp) const;

  void ProjectParticleOnMAC(int idx);

  varType dhat(varType r) const;
private:

};