#include "PIC.hpp"
#include <random>

void PIC::RefillParticles() {
  const int targetPPC = params.ppcx * params.ppcy;

  // single-threaded per cell each Add goes directly into the right cell,
  // no lock needed since we are not parallel here
  for (int cj = 0; cj < ny; cj++) {
    for (int ci = 0; ci < nx; ci++) {
      if (IS_SOLID(fields->Label(ci + 1, cj + 1)))
        continue;

      bool isInflowU = IS_BC_U(fields->Label(ci, cj + 1));
      bool isInflowV = IS_BC_V(fields->Label(ci + 1, cj));
      bool isInflow = isInflowU || isInflowV;

      if (isInflow) {
        // physically coherent emission rate CH7 p.114
        const varType normalVelocity =
            isInflowU ? fields->u.Get(ci, cj) : fields->v.Get(ci, cj);
        const varType speed = std::max(normalVelocity, varType(0));
        if (speed <= REAL_EPSILON)
          continue;

        const varType normalCellWidth = isInflowU ? dx : dy;
        varType dWdt = targetPPC * speed / normalCellWidth;
        varType n = dt * dWdt;

        int toEmit = static_cast<int>(n);
        varType remainder = n - varType(toEmit);
        if (rand01() < remainder)
          toEmit++;

        for (int m = 0; m < toEmit; m++) {
          // Emit on the inflow face. The random birth time below then spreads
          // particles through the swept slab for this step; sampling inside the
          // whole cell as well would double-count advection at the inlet.
          varType x = isInflowU ? ci * dx : (ci + rand01()) * dx;
          varType y = isInflowV ? cj * dy : (cj + rand01()) * dy;

          varType u =
              isInflowU ? fields->u.Get(ci, cj) : fields->interpolateU(x, y);
          varType v =
              isInflowV ? fields->v.Get(ci, cj) : fields->interpolateV(x, y);

          // random birth time + partial advection CH7 p.115
          varType tau = rand01() * dt;
          varType remaining = dt - tau;

          const varType xMax = dx * nx;
          const varType yMax = dy * ny;

          varType xa = std::clamp(x + remaining * u, varType(0),
                                  std::nextafter(xMax, varType(0)));
          varType ya = std::clamp(y + remaining * v, varType(0),
                                  std::nextafter(yMax, varType(0)));

          varType ua = fields->interpolateU(xa, ya);
          varType va = fields->interpolateV(xa, ya);

          const int fi =
              std::clamp(static_cast<int>(std::floor(xa / dx)), 0, nx - 1);
          const int fj =
              std::clamp(static_cast<int>(std::floor(ya / dy)), 0, ny - 1);
          if (IS_SOLID(fields->Label(fi + 1, fj + 1)))
            continue;

          // add directly into the destination cell
          (*cloud)(fi, fj).Add(xa, ya, ua, va, 0);
          fields->setFluid(fi + 1, fj + 1);
        }

      } else if (IS_FLUID(fields->Label(ci + 1, cj + 1))) {
        int alive = cloud->countIn(ci, cj);
        int missing = targetPPC - alive;
        if (missing <= 0)
          continue;

        for (int m = 0; m < missing; m++) {
          varType x = (ci + rand01()) * dx;
          varType y = (cj + rand01()) * dy;

          varType u = fields->interpolateU(x, y);
          varType v = fields->interpolateV(x, y);

          // particle stays in this cell so add directly
          (*cloud)(ci, cj).Add(x, y, u, v, 0);
        }
      }
    }
  }
}
