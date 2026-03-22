#include "Particles.hpp"

// inline ? j'ai des problemes en le rajoutant
void Particles::DropOneParticle(int ip, int jp, 
                    varType x, varType y, varType u, varType v, int id) {
    SetX(ip, jp, x);
    SetY(ip, jp, y);
    SetU(ip, jp, u);
    SetV(ip, jp, v);
    SetId(ip, jp, id);
    SetDead(ip, jp, false);
}

void Particles::InitParticleGrid() {
    unsigned id = 0;

    for (int icell = 0; icell < nx; icell++) {
        for (int jcell = 0; jcell < ny; jcell++) {

            for (int a = 0; a < ppcx; a++) {
                for (int b = 0; b < ppcy; b++) {
                    
                    // PAS REUSSI A FAIRE ça

                    /*if (fields->Label(icell, jcell) == fields->SOLID) { 
                        continue;
                    }*/ 

                    int ip = icell * ppcx + a;
                    int jp = jcell * ppcy + b;

                    varType x = (icell + (a + 0.5) / ppcx) * dx;
                    varType y = (jcell + (b + 0.5) / ppcy) * dy;

                    DropOneParticle(ip, jp, x, y, 0.0, 0.0, id);

                    id++;
                }
            }
        }
    }
}

bool ParticleSlots::PopParticleSlot(int& ip, int& jp)
{
    if (A.empty()) {
        ip = -1;
        jp = -1;
        return false;
    }

    ParticleSlot s = A.back();
    A.pop_back();

    ip = s.ip;
    jp = s.jp;
    return true;
}