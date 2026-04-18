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
    // then if no solid and no particle
    fields->setAir(i + 1, j + 1);
  }
}
    OMP_PRAGMA(omp parallel for collapse(2))
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        if (!IS_AIR(fields->Label(i + 1, j + 1)))
          continue;
        // zero all 4 faces of this cell!
        // Only zero a face if the neighbour on the other side is also non-fluid
        // Left face u(i, j): other side is Label(i, j+1)
        if (!IS_FLUID(fields->Label(i, j + 1)))
          fields->u.Set(i, j, static_cast<varType>(0));
        // Right face u(i+1, j): other side is Label(i+2, j+1)
        if (!IS_FLUID(fields->Label(i + 2, j + 1)))
          fields->u.Set(i + 1, j, static_cast<varType>(0));
        // Bottom face v(i, j): other side is Label(i+1, j)
        if (!IS_FLUID(fields->Label(i + 1, j)))
          fields->v.Set(i, j, static_cast<varType>(0));
        // Top face v(i, j+1): other side is Label(i+1, j+2)
        if (!IS_FLUID(fields->Label(i + 1, j + 2)))
          fields->v.Set(i, j + 1, static_cast<varType>(0));
        // the write race is harmless because we write 0 anyway

        /*varType alive = 0.0;
        // TODO : gerer les BC_U // BC_V // BC_P
        // countAliveParticles is of size (nx, ny)
        // Label is of size (nx + 2; ny + 2) so shifted by i + 1, j + 1
        alive = fields->countAliveParticles.Get(i, j);

        if (alive == 0.0){
            //fields->ResetLabel(i + 1, j + 1, Fields2D::AIR);
            fields->setAir(i+1, j+1);
            // verifier les indices ici dessous
            if(IS_AIR(fields->Label(i, j + 1)))
                fields->u.Set(i, j, 0.0);
            if(IS_AIR(fields->Label(i + 1, j)))
                fields->v.Set(i, j, 0.0);
        }
        else{
            //fields->ResetLabel(i + 1, j + 1, Fields2D::FLUID);
            fields->setFluid(i+1,j+1);
        }
        */
      }
    }
}
