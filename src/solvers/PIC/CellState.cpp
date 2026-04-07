#include "PIC.hpp"

void PIC::UpdateCellState(){
    varType alive = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            // TODO : gerer les BC_U // BC_V // BC_P
            if (fields->Label(i + 1, j + 1) & Fields2D::SOLID) 
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
