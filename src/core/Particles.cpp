#include "Particles.hpp"

void Particles::InitParticleGrid() {
    int id = 0;

    for (int icell = 0; icell < nx; icell++) {
        for (int jcell = 0; jcell < ny; jcell++) {

            for (int a = 0; a < ppcx; a++) {
                for (int b = 0; b < ppcy; b++) {

                    int ip = icell * ppcx + a;
                    int jp = jcell * ppcy + b;

                    varType x = (icell + (a + 0.5) / ppcx) * dx;
                    varType y = (jcell + (b + 0.5) / ppcy) * dy;

                    SetX(ip, jp, x);
                    SetY(ip, jp, y);
                    SetU(ip, jp, 0.0);
                    SetV(ip, jp, 0.0);
                    SetId(ip, jp, id);

                    id++;
                }
            }
        }
    }
}