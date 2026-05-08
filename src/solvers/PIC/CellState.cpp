#include "PIC.hpp"
#include <iostream>

// we can look at the cloud to see if there is particle
// it there is we can simply set fluid
void PIC::UpdateCellState() const {
    OMP_PRAGMA(omp parallel for collapse(2))
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        // early exit
        const labeltype current = fields->Label(i + 1, j + 1);
        if (IS_SOLID(current) || IS_BC_U(current) || IS_BC_V(current))
          continue;
        if (cloud->countIn(i, j) > 0) {
          fields->setFluid(i + 1, j + 1);
          continue;
        }
        fields->setAir(i + 1, j + 1);
      }
    }
    OMP_PRAGMA(omp parallel for collapse(2))
    for (int i = 1; i < fields->p.nx - 1; ++i) {
      for (int j = 1; j < fields->p.ny - 1; ++j) {
        if (!IS_AIR(fields->Label(i, j)))
          continue;

        fields->p.Set(i, j, 0.0);

        labeltype leftLabel = fields->Label(i - 1, j);
        labeltype rightLabel = fields->Label(i + 1, j);
        labeltype bottomLabel = fields->Label(i, j - 1);
        labeltype topLabel = fields->Label(i, j + 1);

        // u, v = 0 at AIR|AIR or AIR|SOLID interfaces
        if (!IS_FLUID(leftLabel))
          fields->u.Set(i - 1, j - 1, static_cast<varType>(0));
        if (!IS_FLUID(rightLabel))
          fields->u.Set(i, j - 1, static_cast<varType>(0));
        if (!IS_FLUID(bottomLabel))
          fields->v.Set(i - 1, j - 1, static_cast<varType>(0));
        if (!IS_FLUID(topLabel))
          fields->v.Set(i - 1, j, static_cast<varType>(0));
      }
    }
}
