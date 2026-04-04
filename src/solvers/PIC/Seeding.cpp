#include "PIC.hpp"
#include <random>

void PIC::RefillParticles() {
  CountAliveParticles();

  const int targetPPC = particles->ppcx * particles->ppcy;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      if (fields->Label(i, j) & Fields2D::SOLID)
        continue;

      int missing =
          targetPPC - static_cast<int>(fields->countAliveParticles.Get(i, j));
      if (missing <= 0)
        continue;

      for (int m = 0; m < missing; m++) {
        int idx = deadSlots->pop();
        if (idx < 0)
          break; // no free slots left this step

        varType x = (i + rand01()) * dx;
        varType y = (j + rand01()) * dy;
        varType u = 0.0, v = 0.0;

        if (fields->Label(i + 1, j) & Fields2D::SOLID) {
              u = 0.0;
              v = interpolateV(x, y);
        } else if (fields->Label(i - 1, j) & Fields2D::SOLID) {
              u = interpolateU(x, y);
              v = 0.0;
        } else {
          u = (fields->Label(i, j) & Fields2D::BC_U) ? 
                    fields->u.Get(i, j) : interpolateU(x, y);

          v = (fields->Label(i, j) & Fields2D::BC_V) ? 
                    fields->v.Get(i, j) : interpolateV(x, y);
        }
        particles->DropOneParticle(idx, x, y, u, v, static_cast<unsigned>(idx));
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
      fields->countAliveParticles.Set(i, j, 0.0);

  const int cap = particles->capacity;
  for (int idx = 0; idx < cap; ++idx) {
    if (particles->IsDead(idx))
      continue;

    int ci = std::clamp(static_cast<int>(std::floor(particles->GetX(idx) / dx)),
                        0, nx - 1);
    int cj = std::clamp(static_cast<int>(std::floor(particles->GetY(idx) / dy)),
                        0, ny - 1);

    if (!(fields->Label(ci, cj) & Fields2D::SOLID)) {
      varType cnt = fields->countAliveParticles.Get(ci, cj);
      fields->countAliveParticles.Set(ci, cj, cnt + 1.0);
    }
  }
}
