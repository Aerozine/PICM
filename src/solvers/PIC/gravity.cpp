#include "PIC.hpp"
/*    must be applied to fields
void PIC::ApplyGravity() const {
  const varType g = params.gravity;
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < fields->v.ny; j++) {
    for (int i = 0; i < fields->v.nx; i++) {
      // v(i,j): bottom cell = Label(i+1, j), top cell = Label(i+1, j+1)
      labeltype bottom = fields->Label(i + 1, j);
      labeltype top = fields->Label(i + 1, j + 1);

      if (IS_SOLID(bottom) || IS_SOLID(top))
        continue;
      if (IS_BC_V(bottom))
        continue;
      // Apply to any fluid-adjacent face, including free-surface faces
      // (original code required both FLUID, starving surface faces of gravity)
      fields->v.Set(i, j, fields->v.Get(i, j) - dt * g);
    }
  }
}
// */

// gravity through particles -> change order of functions in Step() in PIC.cpp
//*
void PIC::ApplyGravity() const {
varType g = 9.81;
OMP_PRAGMA(omp parallel for)
  for (int idx = 0; idx < particles->size(); idx++){
        varType v = particles->GetV(idx);
        particles->SetV(idx, v - dt * g);
  }
}
// */