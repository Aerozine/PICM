#include "PIC.hpp"
#include <cmath>

void PIC::Advect() {
  const varType xMax = dx * nx;
  const varType yMax = dy * ny;
  const int np = particles->size();

  std::vector<uint8_t> keep(np, 1);

  fields->countAliveParticles->reset();

  OMP_PRAGMA(omp parallel for)
  for (int idx = 0; idx < np; ++idx) {
    varType x0 = particles->GetX(idx);
    varType y0 = particles->GetY(idx);
    varType u0 = particles->GetU(idx);
    varType v0 = particles->GetV(idx);

    varType xmid = x0 + static_cast<varType>(0.5) * dt * u0;
    varType ymid = y0 + static_cast<varType>(0.5) * dt * v0;

    varType umid = fields->u.interpolate(xmid, ymid, dx, dy, 0);
    varType vmid = fields->v.interpolate(xmid, ymid, dx, dy, 1);

    varType x1 = x0 + dt * umid;
    varType y1 = y0 + dt * vmid;

    if (x1 < static_cast<varType>(0) || x1 >= xMax ||
        y1 < static_cast<varType>(0) || y1 >= yMax) {
      keep[idx] = 0;
      continue;
    }

    int i1 = std::clamp(static_cast<int>(std::floor(x1 / dx)), 0, nx - 1);
    int j1 = std::clamp(static_cast<int>(std::floor(y1 / dy)), 0, ny - 1);

    if (IS_SOLID(fields->Label(i1 + 1, j1 + 1))) {
      continue; // particle stays at old position, keep[idx] remains 1
    }

    particles->SetX(idx, x1);
    particles->SetY(idx, y1);
  }

  // Single sequential pass: remove dead particles + count alive ones
  for (int idx = np - 1; idx >= 0; --idx) {
    if (!keep[idx]) {
      particles->Remove(idx);
    } else {
      int ci = std::clamp(static_cast<int>(std::floor(particles->GetX(idx) / dx)), 0, nx - 1);
      int cj = std::clamp(static_cast<int>(std::floor(particles->GetY(idx) / dy)), 0, ny - 1);
      if (!IS_SOLID(fields->Label(ci + 1, cj + 1))) {
        varType cnt = fields->countAliveParticles->Get(ci, cj);
        fields->countAliveParticles->Set(ci, cj, cnt + 1.0);
      }
    }
  }
}
