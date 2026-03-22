#include "PIC.hpp"
#include <random>

void PIC::RefillParticles() {
    CountAliveParticles();

    const int targetPPC = 4;
    int ip, jp, id;
    varType x, y, u, v;

    // refill underpopulated fluid cells
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {

            if (fields->Label(i,j) &! Fields2D::FLUID) continue;

            int missing = targetPPC - fields->countAliveParticles.Get(i,j);
            if (missing <= 0) continue;

            for (int m = 0; m < missing; m++) {
                
                if (deadSlots->Empty()) { // injected as many particles as dead ones
                    return; // TODO reinject more ?
                } else if (deadSlots->PopParticleSlot(ip, jp))
                {
                    x = (i + rand01()) * dx;
                    y = (j + rand01()) * dy;

                    if (fields->Label(i, j) & Fields2D::BC_U){
                        u = fields->u.Get(i, j); // !! pas sur concordance indices
                    } else {
                        u = interpolateU(x, y);
                    }

                    if (fields->Label(i, j) & Fields2D::BC_V){
                        v = fields->u.Get(i, j); // !! pas sur concordance indices
                    } else {
                        v = interpolateV(x, y);
                    }

                    id = particles->GetId(ip, jp);

                    particles->DropOneParticle(ip, jp, x, y, u, v, id);
                } else { return; }         
            }
        }
    }
}

varType PIC::rand01() {

    static std::mt19937 rng(std::random_device{}());
    static std::uniform_real_distribution<varType> dist(varType(0), varType(1));
    return dist(rng);
}


void PIC::CountAliveParticles() {

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            fields->countAliveParticles.Set(i,j, 0.0);

    varType counter = 0.0;

    for (int icell = 0; icell < nx; icell++) {
        for (int jcell = 0; jcell < ny; jcell++) {
            for (int a = 0; a < particles->ppcx; ++a) {
                for (int b = 0; b < particles->ppcy; ++b) {
                    int ip = icell * particles->ppcx + a;
                    int jp = jcell * particles->ppcy + b;

                    if (particles->IsDead(ip, jp)) continue;

                    varType x = particles->GetX(ip, jp);
                    varType y = particles->GetY(ip, jp);

                    int i = std::max(0, std::min(nx - 1, int(std::floor(x / dx))));
                    int j = std::max(0, std::min(ny - 1, int(std::floor(y / dy))));

                    if (fields->Label(i,j) == Fields2D::FLUID) {
                        counter = fields->countAliveParticles.Get(i, j); 
                        fields->countAliveParticles.Set(i, j, counter + 1.0);
                    }
                }
            }
        }
    }
}
 