#include "PIC.hpp"

void PIC::ApplyGravity() {
OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < fields->v.ny; j++){
    for (int i = 0; i < fields->v.nx; i++){
      // IF solid underneath -> do not apply gravity
      // i.e. impermeability
      if (IS_SOLID(fields->Label(i + 1, j)))
        continue;
      if (IS_FLUID(fields->Label(i + 1, j + 1)) ||
          IS_FLUID(fields->Label(i + 1, j))){
        varType g = 9.81;
        varType v = fields->v.Get(i ,j);
        fields->v.Set(i ,j, v - dt * g);
      }
    }
  }
}