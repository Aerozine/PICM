#include "PIC.hpp"
#include <cmath>
// @todo can be merged with SemiLagrangian.cpp similarly to Project.cpp
void PIC::AdvectParticles() {
  const varType xMax = dx * nx;
  const varType yMax = dy * ny;
  const int np = particles->size();

  // std::vector<char> keep(np, 1);

OMP_PRAGMA(omp parallel for)
for (int idx = 0; idx < np; ++idx) {
  varType x0 = particles->GetX(idx);
  varType y0 = particles->GetY(idx);
  varType u0 = particles->GetU(idx);
  varType v0 = particles->GetV(idx);

  // varType xmid = x0 + varType(0.5) * dt * u0;
  // varType ymid = y0 + varType(0.5) * dt * v0;
  varType xmid = std::clamp(x0 + varType(0.5) * dt * u0, varType(0),
                            std::nextafter(xMax, varType(0)));
  varType ymid = std::clamp(y0 + varType(0.5) * dt * v0, varType(0),
                            std::nextafter(yMax, varType(0)));
  //@todo CFL and adapt timestep
  varType umid = interpolateU(xmid, ymid);
  varType vmid = interpolateV(xmid, ymid);

  varType x1 = x0 + dt * umid;
  varType y1 = y0 + dt * vmid;
  /*
   *@todo why ?
   *@todo assert particles that goes in with assert
   if (x1 < varType(0) || x1 >= xMax || y1 < varType(0) || y1 >= yMax) {
     keep[idx] = 0;
     continue;
   }
   */
  x1 = std::clamp(x1, varType(0), std::nextafter(xMax, varType(0)));
  y1 = std::clamp(y1, varType(0), std::nextafter(yMax, varType(0)));

  int i1 = std::clamp(static_cast<int>(std::floor(x1 / dx)), 0, nx - 1);
  int j1 = std::clamp(static_cast<int>(std::floor(y1 / dy)), 0, ny - 1);

  // if (fields->Label(i1 + 1, j1 + 1) & Fields2D::SOLID) {
  //  dont move !
  if (IS_SOLID(fields->Label(i1 + 1, j1 + 1))) {
    // keep[idx] = 0;
    // keep particles
    continue;
  }

  particles->SetX(idx, x1);
  particles->SetY(idx, y1);
}
/*
for (int idx = np - 1; idx >= 0; --idx) {
  if (!keep[idx]) {
    particles->Remove(idx);
  }

}
*/
}
