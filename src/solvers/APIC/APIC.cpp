#include "APIC.hpp"

#include <algorithm>
#include <cmath>

APIC::APIC(Parameters &params) : PIC(params) {}

void APIC::ProjectParticlesOnGrid() {
  const int radius = params.kernelOrder;

  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
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

      const varType xFace = static_cast<varType>(i) * dx;
      const varType yFace =
          (static_cast<varType>(j) + varType(0.5)) * dy;

      varType sum = 0;
      varType wt = 0;

      const int ci_lo = std::max(0, i - radius);
      const int ci_hi = std::min(nx - 1, i + radius - 1);
      const int cj_lo = std::max(0, j - radius + 1);
      const int cj_hi = std::min(ny - 1, j + radius);

      for (int ci = ci_lo; ci <= ci_hi; ++ci) {
        for (int cj = cj_lo; cj <= cj_hi; ++cj) {
          const Particles &cell = (*cloud)(ci, cj);
          for (int p = 0; p < cell.size(); ++p) {
            const varType x = cell.GetX(p);
            const varType y = cell.GetY(p);
            const varType xg = x / dx;
            const varType yg = y / dy - varType(0.5);
            const varType k = hat(xg - varType(i)) * hat(yg - varType(j));
            if (k <= varType(0))
              continue;

            const varType ox = xFace - x;
            const varType oy = yFace - y;
            const varType affine =
                cell.GetU(p) + cell.GetCuX(p) * ox + cell.GetCuY(p) * oy;

            sum += k * affine;
            wt += k;
          }
        }
      }

      fields->u.Set(i, j, wt >= REAL_EPSILON ? sum / wt : varType(0));
    }
  }

  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
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

      const varType xFace =
          (static_cast<varType>(i) + varType(0.5)) * dx;
      const varType yFace = static_cast<varType>(j) * dy;

      varType sum = 0;
      varType wt = 0;

      const int ci_lo = std::max(0, i - radius + 1);
      const int ci_hi = std::min(nx - 1, i + radius);
      const int cj_lo = std::max(0, j - radius);
      const int cj_hi = std::min(ny - 1, j + radius - 1);

      for (int ci = ci_lo; ci <= ci_hi; ++ci) {
        for (int cj = cj_lo; cj <= cj_hi; ++cj) {
          const Particles &cell = (*cloud)(ci, cj);
          for (int p = 0; p < cell.size(); ++p) {
            const varType x = cell.GetX(p);
            const varType y = cell.GetY(p);
            const varType xg = x / dx - varType(0.5);
            const varType yg = y / dy;
            const varType k = hat(xg - varType(i)) * hat(yg - varType(j));
            if (k <= varType(0))
              continue;

            const varType ox = xFace - x;
            const varType oy = yFace - y;
            const varType affine =
                cell.GetV(p) + cell.GetCvX(p) * ox + cell.GetCvY(p) * oy;

            sum += k * affine;
            wt += k;
          }
        }
      }

      fields->v.Set(i, j, wt >= REAL_EPSILON ? sum / wt : varType(0));
    }
  }
}

void APIC::accumulateAffineComponent(const Grid2D &grid, varType xg, varType yg,
                                     int imax, int jmax, varType &value,
                                     varType &gradX, varType &gradY) const {
  const int radius = params.kernelOrder;
  const int i0 = static_cast<int>(std::floor(xg));
  const int j0 = static_cast<int>(std::floor(yg));

  value = 0;
  gradX = 0;
  gradY = 0;

  for (int dj = -radius; dj <= radius; ++dj) {
    for (int di = -radius; di <= radius; ++di) {
      const int i = i0 + di;
      const int j = j0 + dj;
      if (i < 0 || i >= imax || j < 0 || j >= jmax)
        continue;

      const varType wx = hat(xg - varType(i));
      const varType wy = hat(yg - varType(j));
      const varType w = wx * wy;
      if (w <= varType(0))
        continue;

      const varType sample = grid.Get(i, j);
      value += w * sample;
      gradX += (dhat(xg - varType(i)) / dx) * wy * sample;
      gradY += wx * (dhat(yg - varType(j)) / dy) * sample;
    }
  }
}

void APIC::ProjectGridOnParticles() {
  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int ci = 0; ci < nx; ++ci) {
    for (int cj = 0; cj < ny; ++cj) {
      Particles &cell = (*cloud)(ci, cj);
      for (int p = 0; p < cell.size(); ++p) {
        const varType x = cell.GetX(p);
        const varType y = cell.GetY(p);

        varType u = 0;
        varType cuX = 0;
        varType cuY = 0;
        accumulateAffineComponent(fields->u, x / dx, y / dy - varType(0.5),
                                  fields->u.nx, fields->u.ny, u, cuX, cuY);

        varType v = 0;
        varType cvX = 0;
        varType cvY = 0;
        accumulateAffineComponent(fields->v, x / dx - varType(0.5), y / dy,
                                  fields->v.nx, fields->v.ny, v, cvX, cvY);

        cell.SetU(p, u);
        cell.SetV(p, v - dt * params.gravity);
        cell.SetCu(p, cuX, cuY);
        cell.SetCv(p, cvX, cvY);
      }
    }
  }
}
