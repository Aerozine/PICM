#include "PIC.hpp"

void PIC::UpdateCellState() const {
  // IN countalive particle cell referential
OMP_PRAGMA(omp parallel for collapse(2))
for (int i = 0; i < fields->countAliveParticles->nx; i++) {
  for (int j = 0; j < fields->countAliveParticles->ny; j++) {
    // early exit
    const labeltype current = fields->Label(i + 1, j + 1);
    if (IS_SOLID(current) || IS_BC_U(current) || IS_BC_V(current))
      continue;
    if (fields->countAliveParticles->Get(i, j) > 0.0) {
      fields->setFluid(i + 1, j + 1);
      continue;
    }
    fields->setAir(i + 1, j + 1);
  }
}
    OMP_PRAGMA(omp parallel for collapse(2))
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        if (!IS_AIR(fields->Label(i + 1, j + 1)))
          continue;
       
        labeltype leftLabel = fields->Label(i, j + 1);
        labeltype rightLabel = fields->Label(i + 2, j + 1);
        labeltype bottomLabel = fields->Label(i + 1, j);
        labeltype topLabel = fields->Label(i + 1, j + 2);

        // u, v = 0 at AIR|AIR or AIR|SOLID interfaces
        if (!IS_FLUID(leftLabel))
          fields->u.Set(i, j, static_cast<varType>(0));
        if (!IS_FLUID(rightLabel))
          fields->u.Set(i + 1, j, static_cast<varType>(0));
        if (!IS_FLUID(bottomLabel))
          fields->v.Set(i, j, static_cast<varType>(0));
        if (!IS_FLUID(topLabel))
          fields->v.Set(i, j + 1, static_cast<varType>(0));
        
        // if (IS_FLUID(leftLabel) &&  IS_FLUID(rightLabel) &&
        //     IS_FLUID(bottomLabel) && IS_FLUID(topLabel))
        //   fields->setFluid(i + 1, j + 1);
      }
    }
}
