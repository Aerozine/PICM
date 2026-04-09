#include "PIC.hpp"

void PIC::ApplyGravity() {
  for (int j = 0; j < fields->v.ny; j++){
    for (int i = 0; i < fields->v.nx; i++){
      if (IS_FLUID(fields->Label(i + 1, j + 1))){
        varType g = 9.81;
        varType v = fields->v.Get(i ,j);
        fields->v.Set(i ,j, v - dt * g);
      }
    }
  }
}