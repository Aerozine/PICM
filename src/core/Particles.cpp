#include "Particles.hpp"

void Particles::DropOneParticle(int idx, varType x, varType y, varType u,
                                varType v, unsigned id) {
  A[idx].pos.x = x;
  A[idx].pos.y = y;
  A[idx].vel.x = u;
  A[idx].vel.y = v;
  A[idx].id = id;
  A[idx].dead = false;
}

void Particles::InitParticleGrid() {
  // Activate the first nx*ny*ppcx*ppcy slots and place them on the
  // regular sub-cell grid.  Remaining slots stay dead (default-constructed).
  unsigned id = 0;
  int idx = 0;

  for (int icell = 0; icell < nx; icell++) {
    for (int jcell = 0; jcell < ny; jcell++) {
      for (int a = 0; a < ppcx; a++) {
        for (int b = 0; b < ppcy; b++) {
          varType x = (icell + (a + varType(0.5)) / ppcx) * dx;
          varType y = (jcell + (b + varType(0.5)) / ppcy) * dy;
          DropOneParticle(idx, x, y, 0.0, 0.0, id);
          ++idx;
          ++id;
        }
      }
    }
  }
  // Slots idx .. capacity-1 remain dead.
}
