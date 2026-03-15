#include "PIC.hpp"

void PIC::RefillParticles() {
    unsigned numSeeds = 2 * dy * ny; // TODO : handle better
    varType spacing = dy * ny / numSeeds;
    
    unsigned seeds = 0;

    varType x, y;
    int id;

    int ppcx = particles->ppcx;
    int ppcy = particles->ppcy;

    for (int icell = 0; icell < nx; icell++) {
        for (int jcell = 0; jcell < ny; jcell++) {

            for (int a = 0; a < ppcx; a++) {
                for (int b = 0; b < ppcy; b++) {
                    int ip = icell * ppcx + a;
                    int jp = jcell * ppcy + b;

                    if (particles->IsDead(ip, jp)) {
                        id = particles->GetId(ip, jp);
                        x = 0.0;
                        y = spacing * seeds;

                        particles->DropOneParticle(ip, jp, x, y, 0.0, 0.0, id);
                        seeds ++;
                        if (seeds >= numSeeds) {
                            return;
                        }
                    }
                }
            }
        }
    }
}  
   
 