#include "PIC.hpp"
#include <iostream>
void PIC::ProjectParticlesOnGrid() {

#ifdef USE_SPEED
  constexpr int R = 1;
#else
  constexpr int R = 2;
#endif

    OMP_PRAGMA(omp parallel for collapse(2) schedule(dynamic,8) )
    for (int j = 0; j < fields->u.ny; ++j) {
      for (int i = 0; i < fields->u.nx; ++i) {

        const labeltype left = fields->Label(i, j + 1);
        const labeltype right = fields->Label(i + 1, j + 1);
        if (IS_SOLID(left) || IS_SOLID(right)) {
          fields->u.Set(i, j, varType(0));
          continue;
        }
        if (IS_BC_U(left))
          continue;

        varType sum = 0, wt = 0;

#ifdef USE_SPEED
        const int ci_lo = std::max(0, i - 1);
        const int ci_hi = std::min(nx - 1, i);
        const int cj_lo = std::max(0, j - 1);
        const int cj_hi = std::min(ny - 1, j + 1);
#else
        const int ci_lo = std::max(0, i - R);
        const int ci_hi = std::min(nx - 1, i + R - 1);
        const int cj_lo = std::max(0, j - R + 1);
        const int cj_hi = std::min(ny - 1, j + R);
#endif

        for (int ci = ci_lo; ci <= ci_hi; ++ci) {
          for (int cj = cj_lo; cj <= cj_hi; ++cj) {
            const Particles &cell = (*cloud)(ci, cj);
            const int particleCount = cell.size();
            for (int p = 0; p < particleCount; ++p) {
              const varType xg = cell.GetX(p) / dx;
              const varType yg = cell.GetY(p) / dy - varType(0.5);
              const varType k = hat(xg - varType(i)) * hat(yg - varType(j));
              if (k > varType(0)) {
                sum += k * cell.GetU(p);
                wt += k;
              }
            }
          }
        }

        fields->u.Set(i, j, wt >= REAL_EPSILON ? sum / wt : varType(0));
      }
    }
    OMP_PRAGMA(omp parallel for collapse(2) schedule(dynamic,8))
    for (int j = 0; j < fields->v.ny; ++j) {
      for (int i = 0; i < fields->v.nx; ++i) {

        const labeltype bottom = fields->Label(i + 1, j);
        const labeltype top = fields->Label(i + 1, j + 1);
        if (IS_SOLID(bottom) || IS_SOLID(top)) {
          fields->v.Set(i, j, varType(0));
          continue;
        }
        if (IS_BC_V(bottom))
          continue;

        varType sum = 0, wt = 0;

#ifdef USE_SPEED
        const int ci_lo = std::max(0, i - 1);
        const int ci_hi = std::min(nx - 1, i + 1);
        const int cj_lo = std::max(0, j - 1);
        const int cj_hi = std::min(ny - 1, j);
#else
        const int ci_lo = std::max(0, i - R + 1);
        const int ci_hi = std::min(nx - 1, i + R);
        const int cj_lo = std::max(0, j - R);
        const int cj_hi = std::min(ny - 1, j + R - 1);
#endif

        for (int ci = ci_lo; ci <= ci_hi; ++ci) {
          for (int cj = cj_lo; cj <= cj_hi; ++cj) {
            const Particles &cell = (*cloud)(ci, cj);
            const int particleCount = cell.size();
            for (int p = 0; p < particleCount; ++p) {
              const varType xg = cell.GetX(p) / dx - varType(0.5);
              const varType yg = cell.GetY(p) / dy;
              const varType k = hat(xg - varType(i)) * hat(yg - varType(j));
              if (k > varType(0)) {
                sum += k * cell.GetV(p);
                wt += k;
              }
            }
          }
        }

        fields->v.Set(i, j, wt >= REAL_EPSILON ? sum / wt : varType(0));
      }
    }
}

void PIC::ProjectGridOnParticles() {
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int ci = 0; ci < nx; ++ci) {
      for (int cj = 0; cj < ny; ++cj) {
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
