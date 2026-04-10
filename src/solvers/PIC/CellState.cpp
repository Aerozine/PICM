#include "PIC.hpp"

void PIC::UpdateCellState(){
OMP_PRAGMA(omp parallel for collapse(2))
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            varType alive = 0.0;
            // TODO : gerer les BC_U // BC_V // BC_P 
            // countAliveParticles is of size (nx, ny)
            // Label is of size (nx + 2; ny + 2) so shifted by i + 1, j + 1
            if (IS_SOLID(fields->Label(i + 1, j + 1))) 
                continue;
            alive = fields->countAliveParticles.Get(i, j);

            if (alive == 0.0){
                fields->ResetLabel(i + 1, j + 1, Fields2D::AIR);
            }
            else{
                fields->ResetLabel(i + 1, j + 1, Fields2D::FLUID);
            }
        }
    }
}
