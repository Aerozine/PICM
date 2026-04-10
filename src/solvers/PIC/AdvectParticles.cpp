#include "PIC.hpp"
#include <cmath>

void PIC::AdvectParticles() {
  const varType xMax = dx * nx; // not dx * (nx - 1) ?
  const varType yMax = dy * ny;

OMP_PRAGMA(omp parallel for)
  for (int idx = 0; idx < particles->size(); ++idx) {

    varType x0 = particles->GetX(idx);
    varType y0 = particles->GetY(idx);
    // varType u0 = particles->GetU(idx); // should not change anything ?
    // varType v0 = particles->GetV(idx); // however it gives weird results 
    varType u0 = interpolateU(x0, y0);    // (cube in free wall becomes a line)
    varType v0 = interpolateV(x0, y0); 


    // RK2 mid-point
    varType xmid = x0 + varType(0.5) * dt * u0;
    varType ymid = y0 + varType(0.5) * dt * v0;
    varType umid = interpolateU(xmid, ymid);
    varType vmid = interpolateV(xmid, ymid);

    varType x1 = x0 + dt * umid;
    varType y1 = y0 + dt * vmid;

    // Remove particles instead ?
    if (x1 < varType(0) || x1 >= xMax || y1 < varType(0) || y1 >= yMax) {
      x1 = std::clamp(x1, varType(0.001)*dx, xMax - varType(0.001)*dx); 
      y1 = std::clamp(y1, varType(0.001)*dy, yMax - varType(0.001)*dy); 
      particles->SetX(idx, x1);
      particles->SetY(idx, y1);
      continue;
    }

    int i1 = std::clamp(static_cast<int>(std::floor(x1 / dx)), 0, nx - 1);
    int j1 = std::clamp(static_cast<int>(std::floor(y1 / dy)), 0, ny - 1);

    // label is shifted -> go at i + 1, j + 1
    if (fields->Label(i1 + 1, j1 + 1) & Fields2D::SOLID) {
      x1 = std::clamp(x1, varType(0.001)*dx, xMax - varType(0.001)*dx); // Ce que
      y1 = std::clamp(y1, varType(0.001)*dy, yMax - varType(0.001)*dy); // propose l'IA
      particles->SetX(idx, x1);
      particles->SetY(idx, y1);
      continue;
    }

    particles->SetX(idx, x1);
    particles->SetY(idx, y1);
  }
}

