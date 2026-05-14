#include "APIC.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

// FIX 2: Don't skip solid-adjacent faces — they hold valid zero velocity.
// Skipping them asymmetrically drops stencil contributions near walls,
// which corrupts the affine matrix at free-surface/wall junctions.
[[nodiscard]] bool sampleFaceForAffine(const Fields2D &fields,
                                       const Grid2D &grid, int i, int j,
                                       bool uComponent, varType &sample) {
  i = std::clamp(i, 0, grid.nx - 1);
  j = std::clamp(j, 0, grid.ny - 1);

  sample = grid.Get(i, j);
  return true;
}

} // namespace

APIC::APIC(Parameters &params) : PIC(params) {}

void APIC::ProjectParticlesOnGrid() {
  const int radius = params.kernelOrder;

  // --- U faces: interior (same loop as before, correct) ---
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
      const varType yFace = (static_cast<varType>(j) + varType(0.5)) * dy;

      varType sum = 0;
      varType wt = 0;

      const int ci_lo = std::max(0, i - radius);
      const int ci_hi = std::min(nx - 1, i + radius - 1);
      const int cj_lo = std::max(0, j - radius + 1);
      const int cj_hi = std::min(ny - 1, j + radius);

      for (int cj = cj_lo; cj <= cj_hi; ++cj) {
        for (int ci = ci_lo; ci <= ci_hi; ++ci) {
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

  // FIX 3: Mirror PIC's explicit right-boundary pass for u (i = nx).
  // Without this, the rightmost u-column is never written from particles.
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

    const varType xFace = static_cast<varType>(nx) * dx;
    const varType yFace = (static_cast<varType>(cj) + varType(0.5)) * dy;

    varType sum = 0;
    varType wt = 0;

    const int ci_lo = std::max(0, nx - radius);
    const int ci_hi = nx - 1;
    const int cj_lo = std::max(0, cj - radius + 1);
    const int cj_hi = std::min(ny - 1, cj + radius);

    for (int nj = cj_lo; nj <= cj_hi; ++nj) {
      for (int ci = ci_lo; ci <= ci_hi; ++ci) {
        const Particles &cell = (*cloud)(ci, nj);
        for (int p = 0; p < cell.size(); ++p) {
          const varType x = cell.GetX(p);
          const varType y = cell.GetY(p);
          const varType xg = x / dx;
          const varType yg = y / dy - varType(0.5);
          const varType k = hat(xg - varType(nx)) * hat(yg - varType(cj));
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

    fields->u.Set(nx, cj, wt >= REAL_EPSILON ? sum / wt : varType(0));
  }

  // --- V faces: interior ---
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

      const varType xFace = (static_cast<varType>(i) + varType(0.5)) * dx;
      const varType yFace = static_cast<varType>(j) * dy;

      varType sum = 0;
      varType wt = 0;

      const int ci_lo = std::max(0, i - radius + 1);
      const int ci_hi = std::min(nx - 1, i + radius);
      const int cj_lo = std::max(0, j - radius);
      const int cj_hi = std::min(ny - 1, j + radius - 1);

      for (int cj = cj_lo; cj <= cj_hi; ++cj) {
        for (int ci = ci_lo; ci <= ci_hi; ++ci) {
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

  // FIX 3: Mirror PIC's explicit top-boundary pass for v (j = ny).
  // Without this, the topmost v-row is never written from particles.
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

    const varType xFace = (static_cast<varType>(ci) + varType(0.5)) * dx;
    const varType yFace = static_cast<varType>(ny) * dy;

    varType sum = 0;
    varType wt = 0;

    const int ci_lo = std::max(0, ci - radius + 1);
    const int ci_hi = std::min(nx - 1, ci + radius);
    const int cj_lo = std::max(0, ny - radius);
    const int cj_hi = ny - 1;

    for (int cj = cj_lo; cj <= cj_hi; ++cj) {
      for (int ni = ci_lo; ni <= ci_hi; ++ni) {
        const Particles &cell = (*cloud)(ni, cj);
        for (int p = 0; p < cell.size(); ++p) {
          const varType x = cell.GetX(p);
          const varType y = cell.GetY(p);
          const varType xg = x / dx - varType(0.5);
          const varType yg = y / dy;
          const varType k = hat(xg - varType(ci)) * hat(yg - varType(ny));
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

    fields->v.Set(ci, ny, wt >= REAL_EPSILON ? sum / wt : varType(0));
  }
}

void APIC::accumulateAffineComponent(const Grid2D &grid, varType xg, varType yg,
                                     int imax, int jmax, bool uComponent,
                                     varType &value, varType &gradX,
                                     varType &gradY) const {
  const int radius = params.kernelOrder;
  const int i0 = static_cast<int>(std::floor(xg));
  const int j0 = static_cast<int>(std::floor(yg));

  varType valueRaw = 0;
  varType gradXRaw = 0;
  varType gradYRaw = 0;
  varType weightSum = 0;

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

      varType sample = 0;
      if (!sampleFaceForAffine(*fields, grid, i, j, uComponent, sample))
        continue;

      // FIX 1: Use unnormalized gradient (standard APIC B-matrix).
      // The quotient-rule correction (-valueRaw * dWeightX / W^2) is
      // physically wrong for the affine matrix: it introduces large spurious
      // gradients near free surfaces where weightSum drops sharply, and
      // those bad cuX/cuY values then blow up divergence in P2G via the
      // affine extrapolation u + cuX*ox + cuY*oy.
      const varType dwx = (dhat(xg - varType(i)) / dx) * wy;
      const varType dwy = wx * (dhat(yg - varType(j)) / dy);

      valueRaw += w * sample;
      gradXRaw += dwx * sample;
      gradYRaw += dwy * sample;
      weightSum += w;
    }
  }

  if (weightSum <= REAL_EPSILON) {
    value = 0;
    gradX = 0;
    gradY = 0;
    return;
  }

  value  = valueRaw  / weightSum;
  gradX  = gradXRaw  / weightSum;
  gradY  = gradYRaw  / weightSum;
}

void APIC::ProjectGridOnParticles() {
  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int cj = 0; cj < ny; ++cj) {
    for (int ci = 0; ci < nx; ++ci) {
      Particles &cell = (*cloud)(ci, cj);
      for (int p = 0; p < cell.size(); ++p) {
        const varType x = cell.GetX(p);
        const varType y = cell.GetY(p);

        varType u = 0;
        varType cuX = 0;
        varType cuY = 0;
        accumulateAffineComponent(fields->u, x / dx, y / dy - varType(0.5),
                                  fields->u.nx, fields->u.ny, true, u, cuX,
                                  cuY);

        varType v = 0;
        varType cvX = 0;
        varType cvY = 0;
        accumulateAffineComponent(fields->v, x / dx - varType(0.5), y / dy,
                                  fields->v.nx, fields->v.ny, false, v, cvX,
                                  cvY);

        cell.SetU(p, u);
        cell.SetV(p, v - dt * params.gravity);
        cell.SetCu(p, cuX, cuY);
        cell.SetCv(p, cvX, cvY);
      }
    }
  }
}

void APIC::VelocityNormFromParticles() const {
  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int cj = 0; cj < ny; ++cj) {
    for (int ci = 0; ci < nx; ++ci) {
      const Particles &cell = (*cloud)(ci, cj);
      const int n = cell.size();
      if (n == 0) {
        fields->normVelocity.Set(ci, cj, varType(0));
        continue;
      }

      varType sumSq = 0;
      for (int p = 0; p < n; ++p) {
        const varType u = cell.GetU(p);
        const varType v = cell.GetV(p);
        sumSq += u * u + v * v;
      }
      fields->normVelocity.Set(
          ci, cj, std::sqrt(sumSq / static_cast<varType>(n)));
    }
  }
}

void APIC::WriteOutput(int step) const {
  if (step % params.sampling_rate != 0)
    return;

  if (params.write_norm_velocity)
    VelocityNormFromParticles();
  if (params.write_vorticity)
    fields->VorticityCenterGrid();

  Solver::WriteOutput(step);

  if (params.write_particles && particlesWriter) {
    const double timeValue = static_cast<double>(step) * static_cast<double>(dt);
    const bool ok = particlesWriter->writeCloud(*cloud, "particles", timeValue);
    if (!ok)
      std::cerr << "[APIC] Warning: failed to write particles at step " << step
                << '\n';
  }
}
