#include "PIC.hpp"
#include <random>

void PIC::RefillParticles() {

  const int targetPPC = particles->ppcx * particles->ppcy;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      if (fields->Label(i + 1, j + 1) & Fields2D::SOLID)
        continue;

      bool isInflow = (fields->Label(i + 1, j + 1) & Fields2D::BC_U) ||
                      (fields->Label(i + 1, j + 1) & Fields2D::BC_V);

      if (isInflow) {
        // physically coherent emission rate CH7 p.114
        varType speed = varType(0);
        if (fields->Label(i + 1, j + 1) & Fields2D::BC_U)
          speed = std::abs(fields->u.Get(i, j));
        else
          speed = std::abs(fields->v.Get(i, j));

        varType dWdt = targetPPC * speed / dx;
        varType n    = dt * dWdt;

        int toEmit = static_cast<int>(n);
        varType remainder = n - varType(toEmit);
        if (rand01() < remainder) toEmit++;

        for (int m = 0; m < toEmit; m++) {
          varType x = (i + rand01()) * dx;
          varType y = (j + rand01()) * dy;

          varType u = (fields->Label(i + 1, j + 1) & Fields2D::BC_U) ?
                          fields->u.Get(i, j) : interpolateU(x, y);
          varType v = (fields->Label(i + 1, j + 1) & Fields2D::BC_V) ?
                          fields->v.Get(i, j) : interpolateV(x, y);

          // random birth time + partial advection CH7 p 115
          varType tau     = rand01() * dt;  
          varType remaining = dt - tau;     

          // Simple Euler for partial advection
          varType xa = x + remaining * u;   
          varType ya = y + remaining * v;  

          varType ua = interpolateU(xa, ya); // really necessary ?
          varType va = interpolateV(xa, ya);

          // In bounds verification
          int fi = std::clamp(static_cast<int>(std::floor(xa / dx)), 0, nx - 1);
          int fj = std::clamp(static_cast<int>(std::floor(ya / dy)), 0, ny - 1);
          if (xa < varType(0) || xa >= dx * nx || ya < varType(0) || ya >= dy * ny)
            continue;
          if (fields->Label(fi + 1, fj + 1) & Fields2D::SOLID)
            continue;

          particles->Add(xa, ya, ua, va, static_cast<unsigned>(particles->size()));
        }

      } else if (! (fields->Label(i + 1, j + 1) & Fields2D::AIR)){
        int missing = targetPPC -
            static_cast<int>(fields->countAliveParticles.Get(i, j));
        if (missing <= 0)
          continue;

        for (int m = 0; m < missing; m++) {
          varType x = (i + rand01()) * dx;
          varType y = (j + rand01()) * dy;
          varType u = 0.0, v = 0.0;

          // if (fields->Label(i + 1, j) & Fields2D::SOLID) {
          //   u = 0.0;
          //   v = interpolateV(x, y);
          // } else if (fields->Label(i - 1, j) & Fields2D::SOLID) {
          //   u = interpolateU(x, y);
          //   v = 0.0;
          // } else {
          u = interpolateU(x, y);
          v = interpolateV(x, y);
          // }

          particles->Add(x, y, u, v, static_cast<unsigned>(particles->size()));
        }
      }
    }
  }
}

