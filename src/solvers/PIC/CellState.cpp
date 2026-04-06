#include "PIC.hpp"

void PIC::UpdateCellState(){
    varType alive = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            if ((fields->Label(i + 1, j + 1) & Fields2D::FLUID) || 
                ! (fields->Label(i + 1, j + 1) & Fields2D::AIR))
                continue;
            alive = fields->countAliveParticles.Get(i, j);

            if (alive == 0.0)
                fields->SetLabel(i + 1, j + 1, Fields2D::AIR);
            else  
                fields->SetLabel(i + 1, j + 1, Fields2D::FLUID);
        }
    }
}

// void PIC::UpdateCellState(){
//     for (int i = 0; i < nx; i++) {
//         for (int j = 0; j < ny; j++) {
//             // ne pas toucher aux solides et aux BC
//             if (fields->Label(i, j) & Fields2D::SOLID)
//                 continue;

//             varType alive = fields->countAliveParticles.Get(i, j);

//             if (alive == 0)
//                 fields->SetLabel(i, j, Fields2D::AIR);
//             else
//                 fields->SetLabel(i, j, Fields2D::FLUID);
//         }
//     }
// }