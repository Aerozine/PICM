#include "Particles.hpp"
#include <iostream>

void Particles::InitParticleGrid(Fields2D &fields) {
    A.clear();
    // A.reserve(nx * ny * ppcx * ppcy);
    unsigned id = 0;

    for (int icell = 0; icell < nx; icell++) {
        for (int jcell = 0; jcell < ny; jcell++) {
            if (fields.Label(icell + 1, jcell + 1) & Fields2D::AIR ||
                fields.Label(icell + 1, jcell + 1) & Fields2D::SOLID)
                continue;
            for (int a = 0; a < ppcx; a++) {
                for (int b = 0; b < ppcy; b++) {
                    varType x = (icell + (a + varType(0.5)) / ppcx) * dx;
                    varType y = (jcell + (b + varType(0.5)) / ppcy) * dy;
                    Add(x, y, 0.0, 0.0, id);
                    ++id;
                }
            }
        }
    }
}
