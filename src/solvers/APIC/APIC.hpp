#pragma once
#include "../PIC/PIC.hpp"

class APIC : public PIC {
public:
  explicit APIC(Parameters &params);

protected:
  void ProjectParticlesOnGrid() override;
  void ProjectGridOnParticles() override;

private:
  void accumulateStandardAffineComponent(
      const Grid2D &grid, varType xg, varType yg, int imax, int jmax,
      bool uComponent, varType particleX, varType particleY, varType &value,
      varType &affineX, varType &affineY) const;


  /*
   * Les quatre valeurs existantes stockent directement C_p:
   *
   *   [ GetCuX(p)  GetCuY(p) ]  ligne u
   *   [ GetCvX(p)  GetCvY(p) ]  ligne v
   *
   * La reconstruction standard des papiers est:
   *
   *   B_p = sum_i w_ip v_i (x_i - x_p)^T;
   *   D_p = sum_i w_ip (x_i - x_p)(x_i - x_p)^T;
   *   C_p = B_p * inverse(D_p);
   */
};
