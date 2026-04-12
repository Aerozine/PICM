#include "PIC.hpp"
#include <random>

void PIC::RefillParticles() {

  const int targetPPC = particles->ppcx * particles->ppcy;

  for (int ci = 0; ci < nx; ci++) {
    for (int cj = 0; cj < ny; cj++) {
      // Never touch solid cells
      // Label of domain cell (ci,cj) is at Label(ci+1, cj+1)
      if (IS_SOLID(fields->Label(ci + 1, cj + 1)))
        continue;

      // BC_U lives on the LEFT ghost cell  : Label(ci,   cj+1)
      // BC_V lives on the BOTTOM ghost cell: Label(ci+1, cj  )
      // matching SceneObjects::applyVelocityU -> SetLabel(i, j,   BC_U) with
      // u(i, j-1) matching SceneObjects::applyVelocityV -> SetLabel(i, j-1,
      // BC_V) with v(i-1, j-1)
      bool isInflowU = IS_BC_U(fields->Label(ci, cj + 1));
      bool isInflowV = IS_BC_V(fields->Label(ci + 1, cj));
      bool isInflow = isInflowU || isInflowV;

      if (isInflow) {
        // physically coherent emission rate CH7 p.114
        varType speed = varType(0);
        if (isInflowU)
          // left face of domain cell (ci,cj) is u(ci, cj)
          speed = std::abs(fields->u.Get(ci, cj));
        else
          // bottom face of domain cell (ci,cj) is v(ci, cj)
          speed = std::abs(fields->v.Get(ci, cj));

        varType dWdt = targetPPC * speed / dx;
        varType n = dt * dWdt;

        int toEmit = static_cast<int>(n);
        varType remainder = n - varType(toEmit);
        if (rand01() < remainder)
          toEmit++;

        for (int m = 0; m < toEmit; m++) {
          varType x = (ci + rand01()) * dx;
          varType y = (cj + rand01()) * dy;

          // assign inflow velocity from the BC face, interpolate the other
          // component
          varType u = isInflowU ? fields->u.Get(ci, cj) : interpolateU(x, y);
          varType v = isInflowV ? fields->v.Get(ci, cj) : interpolateV(x, y);

          // random birth time + partial advection CH7 p 115
          varType tau = rand01() * dt;
          varType remaining = dt - tau;

          // clamp before interpolating to avoid out-of-bound positions
          const varType xMax = dx * nx;
          const varType yMax = dy * ny;

          // Simple Euler for partial advection
          varType xa = std::clamp(x + remaining * u, varType(0),
                                  std::nextafter(xMax, varType(0)));
          varType ya = std::clamp(y + remaining * v, varType(0),
                                  std::nextafter(yMax, varType(0)));

          varType ua = interpolateU(xa, ya);
          varType va = interpolateV(xa, ya);

          // In bounds verification
          int fi = std::clamp(static_cast<int>(std::floor(xa / dx)), 0, nx - 1);
          int fj = std::clamp(static_cast<int>(std::floor(ya / dy)), 0, ny - 1);
          if (IS_SOLID(fields->Label(fi + 1, fj + 1)))
            continue;

          particles->Add(xa, ya, ua, va,
                         static_cast<unsigned>(particles->size()));
        }

      } else if (IS_FLUID(fields->Label(ci + 1, cj + 1))) {
        // Only refill interior fluid cells — ones where all 4 neighbours are
        // also fluid. Surface cells are naturally sparse and must not be topped
        // up, otherwise the free surface acts as a particle source.
        // Neighbours in label space: left=(ci, cj+1), right=(ci+2, cj+1),
        //                            bottom=(ci+1, cj), top=(ci+1, cj+2)
        bool isSurface = !IS_FLUID(fields->Label(ci, cj + 1)) ||
                         !IS_FLUID(fields->Label(ci + 2, cj + 1)) ||
                         !IS_FLUID(fields->Label(ci + 1, cj)) ||
                         !IS_FLUID(fields->Label(ci + 1, cj + 2));
        if (isSurface)
          continue;

        // Cell was just promoted to fluid this step (stray particle) — skip.
        int alive = static_cast<int>(fields->countAliveParticles.Get(ci, cj));
        if (alive <= 0)
          continue;

        int missing = targetPPC - alive;
        if (missing <= 0)
          continue;

        for (int m = 0; m < missing; m++) {
          varType x = (ci + rand01()) * dx;
          varType y = (cj + rand01()) * dy;

          // if (fields->Label(i + 1, j) & Fields2D::SOLID) {
          //   u = 0.0;
          //   v = interpolateV(x, y);
          // } else if (fields->Label(i - 1, j) & Fields2D::SOLID) {
          //   u = interpolateU(x, y);
          //   v = 0.0;
          // } else {
          varType u = interpolateU(x, y);
          varType v = interpolateV(x, y);
          // }

          particles->Add(x, y, u, v, static_cast<unsigned>(particles->size()));
        }
      }
    }
  }
}
