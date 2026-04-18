#include "Particles.hpp"
#include <iostream>
#include <random>

varType rand01() {
  static std::mt19937 rng(std::random_device{}());
  static std::uniform_real_distribution<varType> dist(static_cast<varType>(0), static_cast<varType>(1));
  return dist(rng);
}

void Particles::InitParticleGrid(const Fields2D &fields) {
  A.clear();
  // A.reserve(nx * ny * ppcx * ppcy);
  unsigned id = 0;

  for (int icell = 0; icell < nx; icell++) {
    for (int jcell = 0; jcell < ny; jcell++) {
      if (IS_FLUID(fields.Label(icell + 1, jcell + 1))) {
        for (int a = 0; a < ppcx; a++) {
          for (int b = 0; b < ppcy; b++) {
            varType x = (icell + rand01()) * dx;
            varType y = (jcell + rand01()) * dy;
            Add(x, y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, id);
            ++id;
          }
        }
      }
    }
  }
}
