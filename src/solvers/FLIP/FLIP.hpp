/*
#pragma once
#include "../PIC/PIC.hpp"

class FLIP : public PIC {
public:
  explicit FLIP(Parameters &params);

  void Step() override;

protected:
  void ProjectGridOnParticles() override;
  Grid2D u_old;
  Grid2D v_old;

  void SaveOldVelocities();
};
*/