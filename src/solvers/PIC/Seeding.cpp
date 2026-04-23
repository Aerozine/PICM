#include "PIC.hpp"
#include <random>

void PIC::RefillParticles() {

  const int targetPPC = particles->ppcx * particles->ppcy;

  for (int ci = 0; ci < nx; ci++) {
    for (int cj = 0; cj < ny; cj++) {
      if (IS_SOLID(fields->Label(ci + 1, cj + 1)))
        continue;

      bool isInflowU = IS_BC_U(fields->Label(ci, cj + 1));
      bool isInflowV = IS_BC_V(fields->Label(ci + 1, cj));
      bool isInflow = isInflowU || isInflowV;

      if (isInflow) {
        // physically coherent emission rate CH7 p.114
        varType speed = varType(0);
        if (isInflowU)
          speed = std::abs(fields->u.Get(ci, cj));
        else
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

          varType u =
              isInflowU ? fields->u.Get(ci, cj) : fields->u.interpolate(x, y, dx, dy, 0);;
          varType v =
              isInflowV ? fields->v.Get(ci, cj) : fields->v.interpolate(x, y, dx, dy, 1);

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

          varType ua =fields->u.interpolate(xa, ya, dx, dy, 0);
          varType va =fields->v.interpolate(xa, ya, dx, dy, 1);

          // In bounds verification
          int fi = std::clamp(static_cast<int>(std::floor(xa / dx)), 0, nx - 1);
          int fj = std::clamp(static_cast<int>(std::floor(ya / dy)), 0, ny - 1);
          if (IS_SOLID(fields->Label(fi + 1, fj + 1)))
            continue;

          particles->Add(xa, ya, ua, va,
                         static_cast<unsigned>(particles->size()));
        }

      } else if (IS_FLUID(fields->Label(ci + 1, cj + 1))) {
        int alive = static_cast<int>(fields->countAliveParticles->Get(ci, cj));
        if (alive <= 0)
          continue;

        int missing = targetPPC - alive;
        if (missing <= 0)
          continue;

        for (int m = 0; m < missing; m++) {
          varType x = (ci + rand01()) * dx;
          varType y = (cj + rand01()) * dy;

          varType u = fields->u.interpolate(x, y, dx, dy, 0);
            //interpolateU(fields->u, x, y);
          varType v = fields->v.interpolate(x, y, dx, dy, 1);
            //interpolateV(fields->v, x, y);

          particles->Add(x, y, u, v, static_cast<unsigned>(particles->size()));
        }
      }
    }
  }
}
