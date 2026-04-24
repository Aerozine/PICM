#pragma once

#include "../PIC/PIC.hpp"

class APIC : public PIC {
public:
  explicit APIC(Parameters &params);

protected:
  void ProjectParticlesOnGrid() override;
  void ProjectGridOnParticles() override;

private:
  void accumulateAffineComponent(const Grid2D &grid, varType xg, varType yg,
                                 int imax, int jmax, varType &value,
                                 varType &gradX, varType &gradY) const;
};
