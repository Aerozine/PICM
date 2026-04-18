#include "PIC.hpp"
#include <cmath>

void PIC::Advect() {
  const varType xMax = dx * nx;
  const varType yMax = dy * ny;
  const int np = particles->size();

  std::vector<char> keep(np, 1);

OMP_PRAGMA(omp parallel for)
  for (int idx = 0; idx < np; ++idx) {
    // On pourrait calculer CountAliveParticles directement ici
    // apres avoir poussé les particules dans leur nouvelle cell
    varType x0 = particles->GetX(idx);
    varType y0 = particles->GetY(idx);
    varType u0 = particles->GetU(idx);
    varType v0 = particles->GetV(idx);

    varType xmid = x0 + static_cast<varType>(0.5) * dt * u0;
    varType ymid = y0 + static_cast<varType>(0.5) * dt * v0;
    varType umid = interpolateU(fields->u, xmid, ymid);
    varType vmid = interpolateV(fields->v, xmid, ymid);

    varType x1 = x0 + dt * umid;
    varType y1 = y0 + dt * vmid;

    if (x1 < static_cast<varType>(0) || x1 >= xMax || y1 < static_cast<varType>(0) || y1 >= yMax) {
      keep[idx] = 0;
      continue;
    }

    int i1 = std::clamp(static_cast<int>(std::floor(x1 / dx)), 0, nx - 1);
    int j1 = std::clamp(static_cast<int>(std::floor(y1 / dy)), 0, ny - 1);

    // if (IS_SOLID(fields->Label(i1 + 1, j1 + 1))) {
    //   keep[idx] = 0;
    //   continue;
    // }

    if (IS_SOLID(fields->Label(i1 + 1, j1 + 1))) {
      // Rewind along the displacement vector until we find a non-solid cell.
      // Binary search between (x0,y0) [known fluid] and (x1,y1) [solid].
      varType xA = x0, yA = y0;   // safe side
      varType xB = x1, yB = y1;   // solid side

      static constexpr int MAX_BISECT = 16;
      for (int k = 0; k < MAX_BISECT; ++k) {
        const varType xM = static_cast<varType>(0.5) * (xA + xB);
        const varType yM = static_cast<varType>(0.5) * (yA + yB);
        const int iM = std::clamp(static_cast<int>(std::floor(xM / dx)), 0, nx - 1);
        if (int jM = std::clamp(static_cast<int>(std::floor(yM / dy)), 0, ny - 1); IS_SOLID(fields->Label(iM + 1, jM + 1)))
          xB = xM, yB = yM;   // midpoint still solid, shrink toward safe
        else
          xA = xM, yA = yM;   // midpoint is fluid, push safe side forward
      }

      // xA,yA is now within a fluid cell, at most dx/2^MAX_BISECT from the
      // solid interface — well inside one cell for MAX_BISECT >= 6.
      // Zero the normal velocity component to avoid re-penetration next step.
      varType ux = interpolateU(fields->u, xA, yA);
      varType vx = interpolateV(fields->v, xA, yA);

      particles->SetX(idx, xA);
      particles->SetY(idx, yA);
      particles->SetU(idx, ux);
      particles->SetV(idx, vx);
      continue;
    }

    particles->SetX(idx, x1);
    particles->SetY(idx, y1);
  }

  for (int idx = np - 1; idx >= 0; --idx) {
    if (!keep[idx]) {
      particles->Remove(idx);
    }
  }
}

