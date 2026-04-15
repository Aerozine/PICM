#pragma once
#include "PIC.hpp"

class FLIP : public PIC {
public:
  FLIP(const Parameters &params);

  void Step() override;
  void ProjectGridOnParticles() override;

private:
  Grid2D u_old;
  Grid2D v_old;

  void SaveOldVelocities();
};