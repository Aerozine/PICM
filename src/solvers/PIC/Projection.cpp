#include "PIC.hpp"
#include <iostream>
void PIC::ProjectParticlesOnGrid() {

#ifdef USE_SPEED
  constexpr int R = 1;
#else
  constexpr int R = 2;
#endif

  OMP_PRAGMA(omp parallel for collapse(2) )
  for (int cj = 0; cj < ny; ++cj) {
    for (int ci = 0; ci < nx; ++ci) {
      const labeltype uLeft = fields->Label(ci, cj + 1);
      const labeltype uRight = fields->Label(ci + 1, cj + 1);
      const bool keepU = IS_BC_U(uLeft);
      const bool zeroU = !keepU && (IS_SOLID(uLeft) || IS_SOLID(uRight));

      const labeltype vBottom = fields->Label(ci + 1, cj);
      const labeltype vTop = fields->Label(ci + 1, cj + 1);
      const bool keepV = IS_BC_V(vBottom);
      const bool zeroV = !keepV && (IS_SOLID(vBottom) || IS_SOLID(vTop));

      if (zeroU)
        fields->u.Set(ci, cj, varType(0));
      if (zeroV)
        fields->v.Set(ci, cj, varType(0));
      if (keepU && keepV)
        continue;

      const bool gatherU = !keepU && !zeroU;
      const bool gatherV = !keepV && !zeroV;
      if (!gatherU && !gatherV)
        continue;

      const int ni_lo = std::max(0, ci - R);
      const int ni_hi = std::min(nx - 1, ci + R);
      const int nj_lo = std::max(0, cj - R);
      const int nj_hi = std::min(ny - 1, cj + R);

      varType sumU = 0;
      varType weightU = 0;
      varType sumV = 0;
      varType weightV = 0;

      if (gatherU && gatherV) {
        for (int nj = nj_lo; nj <= nj_hi; ++nj) {
          for (int ni = ni_lo; ni <= ni_hi; ++ni) {
            const Particles &cell = (*cloud)(ni, nj);
            const int particleCount = cell.size();
            for (int p = 0; p < particleCount; ++p) {
              const varType xg = cell.GetX(p) / dx;
              const varType yg = cell.GetY(p) / dy;

              const varType ku =
                  hat(xg - varType(ci)) * hat(yg - varType(cj) - varType(0.5));
              if (ku > varType(0)) {
                sumU += ku * cell.GetU(p);
                weightU += ku;
              }

              const varType kv =
                  hat(xg - varType(ci) - varType(0.5)) * hat(yg - varType(cj));
              if (kv > varType(0)) {
                sumV += kv * cell.GetV(p);
                weightV += kv;
              }
            }
          }
        }
      } else if (gatherU) {
        for (int nj = nj_lo; nj <= nj_hi; ++nj) {
          for (int ni = ni_lo; ni <= ni_hi; ++ni) {
            const Particles &cell = (*cloud)(ni, nj);
            const int particleCount = cell.size();
            for (int p = 0; p < particleCount; ++p) {
              const varType xg = cell.GetX(p) / dx;
              const varType yg = cell.GetY(p) / dy;
              const varType ku =
                  hat(xg - varType(ci)) * hat(yg - varType(cj) - varType(0.5));
              if (ku > varType(0)) {
                sumU += ku * cell.GetU(p);
                weightU += ku;
              }
            }
          }
        }
      } else {
        for (int nj = nj_lo; nj <= nj_hi; ++nj) {
          for (int ni = ni_lo; ni <= ni_hi; ++ni) {
            const Particles &cell = (*cloud)(ni, nj);
            const int particleCount = cell.size();
            for (int p = 0; p < particleCount; ++p) {
              const varType xg = cell.GetX(p) / dx;
              const varType yg = cell.GetY(p) / dy;
              const varType kv =
                  hat(xg - varType(ci) - varType(0.5)) * hat(yg - varType(cj));
              if (kv > varType(0)) {
                sumV += kv * cell.GetV(p);
                weightV += kv;
              }
            }
          }
        }
      }

      if (gatherU)
        fields->u.Set(ci, cj,
                      weightU >= REAL_EPSILON ? sumU / weightU : varType(0));
      if (gatherV)
        fields->v.Set(ci, cj,
                      weightV >= REAL_EPSILON ? sumV / weightV : varType(0));
    }
  }

  OMP_PRAGMA(omp parallel for schedule(static))
  for (int cj = 0; cj < ny; ++cj) {
    const labeltype left = fields->Label(nx, cj + 1);
    const labeltype right = fields->Label(nx + 1, cj + 1);
    if (IS_BC_U(left))
      continue;
    if (IS_SOLID(left) || IS_SOLID(right)) {
      fields->u.Set(nx, cj, varType(0));
      continue;
    }

    varType sum = 0;
    varType weight = 0;

#ifdef USE_SPEED
    const int ci_lo = std::max(0, nx - 1);
    const int ci_hi = nx - 1;
    const int cj_lo = std::max(0, cj - 1);
    const int cj_hi = std::min(ny - 1, cj + 1);
#else
    const int ci_lo = std::max(0, nx - R);
    const int ci_hi = nx - 1;
    const int cj_lo = std::max(0, cj - R + 1);
    const int cj_hi = std::min(ny - 1, cj + R);
#endif

    for (int nj = cj_lo; nj <= cj_hi; ++nj) {
      for (int ci = ci_lo; ci <= ci_hi; ++ci) {
        const Particles &cell = (*cloud)(ci, nj);
        const int particleCount = cell.size();
        for (int p = 0; p < particleCount; ++p) {
          const varType xg = cell.GetX(p) / dx;
          const varType yg = cell.GetY(p) / dy - varType(0.5);
          const varType k = hat(xg - varType(nx)) * hat(yg - varType(cj));
          if (k > varType(0)) {
            sum += k * cell.GetU(p);
            weight += k;
          }
        }
      }
    }

    fields->u.Set(nx, cj, weight >= REAL_EPSILON ? sum / weight : varType(0));
  }

  OMP_PRAGMA(omp parallel for schedule(static))
  for (int ci = 0; ci < nx; ++ci) {
    const labeltype bottom = fields->Label(ci + 1, ny);
    const labeltype top = fields->Label(ci + 1, ny + 1);
    if (IS_BC_V(bottom))
      continue;
    if (IS_SOLID(bottom) || IS_SOLID(top)) {
      fields->v.Set(ci, ny, varType(0));
      continue;
    }

    varType sum = 0;
    varType weight = 0;

#ifdef USE_SPEED
    const int ci_lo = std::max(0, ci - 1);
    const int ci_hi = std::min(nx - 1, ci + 1);
    const int cj_lo = std::max(0, ny - 1);
    const int cj_hi = ny - 1;
#else
    const int ci_lo = std::max(0, ci - R + 1);
    const int ci_hi = std::min(nx - 1, ci + R);
    const int cj_lo = std::max(0, ny - R);
    const int cj_hi = ny - 1;
#endif

    for (int cj = cj_lo; cj <= cj_hi; ++cj) {
      for (int ni = ci_lo; ni <= ci_hi; ++ni) {
        const Particles &cell = (*cloud)(ni, cj);
        const int particleCount = cell.size();
        for (int p = 0; p < particleCount; ++p) {
          const varType xg = cell.GetX(p) / dx - varType(0.5);
          const varType yg = cell.GetY(p) / dy;
          const varType k = hat(xg - varType(ci)) * hat(yg - varType(ny));
          if (k > varType(0)) {
            sum += k * cell.GetV(p);
            weight += k;
          }
        }
      }
    }

    fields->v.Set(ci, ny, weight >= REAL_EPSILON ? sum / weight : varType(0));
  }
}

void PIC::ProjectGridOnParticles() {
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int cj = 0; cj < ny; ++cj) {
      for (int ci = 0; ci < nx; ++ci) {
        Particles &cell = (*cloud)(ci, cj);
        for (int p = 0; p < cell.size(); ++p) {
          const varType x = cell.GetX(p);
          const varType y = cell.GetY(p);
          cell.SetU(p, fields->u.interpolate<0>(x, y, dx, dy));
          cell.SetV(p, fields->v.interpolate<1>(x, y, dx, dy) -
                           dt * params.gravity);
        }
      }
    }
}
