#include "PIC.hpp"

// @todo why not simply do -1 when killing a particle ?
// and +1 when creating a new one ?
void PIC::CountAliveParticles() {
OMP_PRAGMA(omp parallel for collapse(2))
for (int i = 0; i < nx; i++)
  for (int j = 0; j < ny; j++)
    fields->countAliveParticles.Set(i, j, 0.0);
for (int idx = 0; idx < particles->size(); ++idx) {
  int ci = std::clamp(static_cast<int>(std::floor(particles->GetX(idx) / dx)),
                      0, nx - 1);
  int cj = std::clamp(static_cast<int>(std::floor(particles->GetY(idx) / dy)),
                      0, ny - 1);

  if (!(fields->Label(ci + 1, cj + 1) & Fields2D::SOLID)) {
    varType cnt = fields->countAliveParticles.Get(ci, cj);
    fields->countAliveParticles.Set(ci, cj, cnt + 1.0);
  }
}
}