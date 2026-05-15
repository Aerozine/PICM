#include "PIC.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

//  Surface tension
//  Label(i,j) : i,j ∈ [0, nx+1] × [0, ny+1]  (ghost cells)
//  phi/kappa/normalX/normalY(i,j) : i,j ∈ [0,nx-1] × [0,ny-1]
//  phi(i,j) <-> Label(i+1, j+1)
//  interface_u(i,j) : (nx+1)×ny - face btw phi(i-1,j) & phi(i,j)
//  interface_v(i,j) : nx×(ny+1) - face btw phi(i,j-1) & phi(i,j)

static constexpr varType CSF_KERNEL_RADIUS_FACTOR = REAL_LITERAL(3.0);
static constexpr int CSF_PHI_SMOOTH_PASSES = 2;
static constexpr varType CSF_NORMAL_GRAD_EPS = REAL_LITERAL(1e-4);
static constexpr varType PI = REAL_LITERAL(3.14159265358979);

static inline varType smoothingKernel(varType r, varType h) {
  const varType q = r / h;
  if (q >= 1.0)
    return 0.0;
  const varType s = 1.0 - q;
  return s * s * s;
}

/* // Muller's smoothing kernel -> reference paper on level sets
static inline varType smoothingKernel(varType r, varType h)
{
    const varType q = r / h;

    if (q >= 1.0) return 0.0;

    const varType s = 1.0 - q * q;

    return s * s * s;
}*/

void PIC::UpdatePhiFromParticles() const {
  const int nx = fields->nx;
  const int ny = fields->ny;
  const varType dx = fields->dx;
  const varType dy = fields->dy;
  const varType cellSize = dx;
  const varType h_smooth = CSF_KERNEL_RADIUS_FACTOR * cellSize;
  const int Rsmooth = static_cast<int>(std::ceil(h_smooth / cellSize)) + 1;
  const int targetPPC = std::max(1, params.ppcx * params.ppcy);
  const varType cellArea =
      std::max(dx * dy, std::numeric_limits<varType>::min());
  const varType particleDensity = static_cast<varType>(targetPPC) / cellArea;
  const varType fullKernelValue =
      std::max(particleDensity * PI * h_smooth * h_smooth / REAL_LITERAL(10.0),
               std::numeric_limits<varType>::min());

  fields->phi->reset();
  fields->kappa->reset();
  fields->normalX->reset();
  fields->normalY->reset();

  std::vector<varType> phiValues(static_cast<std::size_t>(nx) * ny,
                                 REAL_LITERAL(0.5) * h_smooth);
  auto idx = [nx](int i, int j) {
    return static_cast<std::size_t>(nx) * j + i;
  };

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {

        const varType xc = (i + 0.5) * dx;
        const varType yc = (j + 0.5) * dy;

        const int ci_lo = std::max(0, i - Rsmooth);
        const int ci_hi = std::min(nx - 1, i + Rsmooth);
        const int cj_lo = std::max(0, j - Rsmooth);
        const int cj_hi = std::min(ny - 1, j + Rsmooth);

        varType cval = 0.0;
        for (int cj = cj_lo; cj <= cj_hi; ++cj) {
          for (int ci = ci_lo; ci <= ci_hi; ++ci) {
            const Particles &cell = (*cloud)(ci, cj);
            const int np = cell.size();
            for (int p = 0; p < np; ++p) {
              const varType r =
                  std::hypot(cell.GetX(p) - xc, cell.GetY(p) - yc);
              cval += smoothingKernel(r, h_smooth);
            }
          }
        }
        const varType volume = std::clamp(cval / fullKernelValue,
                                          REAL_LITERAL(0.0), REAL_LITERAL(1.0));
        phiValues[idx(i, j)] = (REAL_LITERAL(0.5) - volume) * h_smooth;
      }
    }

    std::vector<varType> smoothed(phiValues.size(), REAL_LITERAL(0.0));
    for (int pass = 0; pass < CSF_PHI_SMOOTH_PASSES; ++pass) {
        OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
        for (int j = 0; j < ny; ++j) {
          for (int i = 0; i < nx; ++i) {
            const labeltype label = fields->Label(i + 1, j + 1);
            if (IS_SOLID(label)) {
              smoothed[idx(i, j)] = phiValues[idx(i, j)];
              continue;
            }

            varType sum = REAL_LITERAL(0.0);
            varType weightSum = REAL_LITERAL(0.0);
            for (int dj = -1; dj <= 1; ++dj) {
              const int jj = std::clamp(j + dj, 0, ny - 1);
              const varType wy =
                  (dj == 0) ? REAL_LITERAL(2.0) : REAL_LITERAL(1.0);
              for (int di = -1; di <= 1; ++di) {
                const int ii = std::clamp(i + di, 0, nx - 1);
                if (IS_SOLID(fields->Label(ii + 1, jj + 1)))
                  continue;
                const varType wx =
                    (di == 0) ? REAL_LITERAL(2.0) : REAL_LITERAL(1.0);
                const varType w = wx * wy;
                sum += w * phiValues[idx(ii, jj)];
                weightSum += w;
              }
            }
            smoothed[idx(i, j)] = weightSum > REAL_EPSILON
                                      ? sum / weightSum
                                      : phiValues[idx(i, j)];
          }
        }
        phiValues.swap(smoothed);
    }

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        fields->phi->Set(i, j, phiValues[idx(i, j)]);
}

void PIC::ComputeSurfaceTensionOnFaces() const {
  const int nx = fields->nx;
  const int ny = fields->ny;
  const varType dx = fields->dx;
  const varType dy = fields->dy;
  const varType cellSize = std::min(dx, dy);
  const varType gamma = params.gamma;
  const varType cosTheta = std::cos(params.contactAngle);
  auto idx = [nx](int i, int j) {
    return static_cast<std::size_t>(nx) * j + i;
  };

  fields->interface_u->reset();
  fields->interface_v->reset();

  auto phiGhostX = [&](int i, int j, int di) -> varType {
    int ni = i + di;
    if (ni < 0 || ni >= nx)
      return fields->phi->Get(i, j);
    if (IS_SOLID(fields->Label(ni + 1, j + 1)))
      return fields->phi->Get(i, j);
    return fields->phi->Get(ni, j);
  };

  auto phiGhostY = [&](int i, int j, int dj) -> varType {
    int nj = j + dj;
    if (nj < 0 || nj >= ny)
      return fields->phi->Get(i, j);
    if (IS_SOLID(fields->Label(i + 1, nj + 1)))
      return fields->phi->Get(i, j);
    return fields->phi->Get(i, nj);
  };

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {

        const varType phi_xp = phiGhostX(i, j, +1);
        const varType phi_xm = phiGhostX(i, j, -1);
        const varType phi_yp = phiGhostY(i, j, +1);
        const varType phi_ym = phiGhostY(i, j, -1);

        const varType gx = (phi_xp - phi_xm) / (2.0 * dx);
        const varType gy = (phi_yp - phi_ym) / (2.0 * dy);
        const varType gn = std::hypot(gx, gy);

        if (gn > CSF_NORMAL_GRAD_EPS) {
          fields->normalX->Set(i, j, gx / gn);
          fields->normalY->Set(i, j, gy / gn);
        }
      }
    }

    const varType targetWallComponent =
        std::clamp(-cosTheta, REAL_LITERAL(-1.0), REAL_LITERAL(1.0));
    const varType targetTangentComponent = std::sqrt(std::max(
        REAL_LITERAL(0.0),
        REAL_LITERAL(1.0) - targetWallComponent * targetWallComponent));

    const auto wallNormal = [&](int i, int j, varType &wx, varType &wy) {
      wx = REAL_LITERAL(0.0);
      wy = REAL_LITERAL(0.0);
      if (i <= 0 || IS_SOLID(fields->Label(i, j + 1)))
        wx += REAL_LITERAL(1.0);
      if (i >= nx - 1 || IS_SOLID(fields->Label(i + 2, j + 1)))
        wx -= REAL_LITERAL(1.0);
      if (j <= 0 || IS_SOLID(fields->Label(i + 1, j)))
        wy += REAL_LITERAL(1.0);
      if (j >= ny - 1 || IS_SOLID(fields->Label(i + 1, j + 2)))
        wy -= REAL_LITERAL(1.0);

      const varType wn = std::hypot(wx, wy);
      if (wn <= REAL_EPSILON)
        return false;
      wx /= wn;
      wy /= wn;
      return true;
    };

    const auto nearFluidAirInterface = [&](int i, int j) {
      bool hasFluid = IS_FLUID(fields->Label(i + 1, j + 1));
      bool hasAir = IS_AIR(fields->Label(i + 1, j + 1));
      for (int dj = -1; dj <= 1; ++dj) {
        for (int di = -1; di <= 1; ++di) {
          const int ii = i + di;
          const int jj = j + dj;
          if (ii < 0 || ii >= nx || jj < 0 || jj >= ny)
            continue;
          const labeltype label = fields->Label(ii + 1, jj + 1);
          if (IS_SOLID(label))
            continue;
          hasFluid = hasFluid || IS_FLUID(label);
          hasAir = hasAir || IS_AIR(label);
        }
      }
      return hasFluid && hasAir;
    };

    const auto airDirection = [&](int i, int j, varType &ax, varType &ay) {
      ax = REAL_LITERAL(0.0);
      ay = REAL_LITERAL(0.0);
      for (int dj = -1; dj <= 1; ++dj) {
        for (int di = -1; di <= 1; ++di) {
          if (di == 0 && dj == 0)
            continue;
          const int ii = i + di;
          const int jj = j + dj;
          if (ii < 0 || ii >= nx || jj < 0 || jj >= ny)
            continue;
          const labeltype label = fields->Label(ii + 1, jj + 1);
          if (IS_AIR(label)) {
            ax += static_cast<varType>(di);
            ay += static_cast<varType>(dj);
          }
        }
      }
    };

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        if (IS_SOLID(fields->Label(i + 1, j + 1)))
          continue;

        varType wx, wy;
        if (!wallNormal(i, j, wx, wy) || !nearFluidAirInterface(i, j))
          continue;

        const varType nx0 = fields->normalX->Get(i, j);
        const varType ny0 = fields->normalY->Get(i, j);
        const varType normalWallComponent = nx0 * wx + ny0 * wy;
        varType tx = nx0 - normalWallComponent * wx;
        varType ty = ny0 - normalWallComponent * wy;
        varType tn = std::hypot(tx, ty);

        if (tn <= CSF_NORMAL_GRAD_EPS) {
          varType ax, ay;
          airDirection(i, j, ax, ay);
          const varType airWallComponent = ax * wx + ay * wy;
          tx = ax - airWallComponent * wx;
          ty = ay - airWallComponent * wy;
          tn = std::hypot(tx, ty);
        }

        if (tn <= CSF_NORMAL_GRAD_EPS) {
          tx = -wy;
          ty = wx;
          tn = REAL_LITERAL(1.0);
        }

        tx /= tn;
        ty /= tn;
        fields->normalX->Set(
            i, j, targetTangentComponent * tx + targetWallComponent * wx);
        fields->normalY->Set(
            i, j, targetTangentComponent * ty + targetWallComponent * wy);
      }
    }

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {

        if (IS_SOLID(fields->Label(i + 1, j + 1)))
          continue;

        auto nxVal = [&](int ii, int jj) -> varType {
          if (ii < 0 || ii >= nx || jj < 0 || jj >= ny)
            return fields->normalX->Get(i, j);
          if (IS_SOLID(fields->Label(ii + 1, jj + 1)))
            return fields->normalX->Get(i, j);
          return fields->normalX->Get(ii, jj);
        };
        auto nyVal = [&](int ii, int jj) -> varType {
          if (ii < 0 || ii >= nx || jj < 0 || jj >= ny)
            return fields->normalY->Get(i, j);
          if (IS_SOLID(fields->Label(ii + 1, jj + 1)))
            return fields->normalY->Get(i, j);
          return fields->normalY->Get(ii, jj);
        };

        const varType dnx_dx = (nxVal(i + 1, j) - nxVal(i - 1, j)) / (2.0 * dx);
        const varType dny_dy = (nyVal(i, j + 1) - nyVal(i, j - 1)) / (2.0 * dy);

        fields->kappa->Set(i, j, -(dnx_dx + dny_dy));
      }
    }

    std::vector<varType> rawKappa(static_cast<std::size_t>(nx) * ny,
                                  REAL_LITERAL(0.0));
    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        rawKappa[idx(i, j)] = fields->kappa->Get(i, j);

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        if (IS_SOLID(fields->Label(i + 1, j + 1)))
          continue;
        varType sum = REAL_LITERAL(0.0);
        varType weightSum = REAL_LITERAL(0.0);
        for (int dj = -1; dj <= 1; ++dj) {
          const int jj = std::clamp(j + dj, 0, ny - 1);
          const varType wy = (dj == 0) ? REAL_LITERAL(2.0) : REAL_LITERAL(1.0);
          for (int di = -1; di <= 1; ++di) {
            const int ii = std::clamp(i + di, 0, nx - 1);
            if (IS_SOLID(fields->Label(ii + 1, jj + 1)))
              continue;
            const varType wx =
                (di == 0) ? REAL_LITERAL(2.0) : REAL_LITERAL(1.0);
            const varType w = wx * wy;
            sum += w * rawKappa[idx(ii, jj)];
            weightSum += w;
          }
        }
        fields->kappa->Set(i, j,
                           weightSum > REAL_EPSILON ? sum / weightSum
                                                    : rawKappa[idx(i, j)]);
      }
    }

    const varType maxAbsKappa =
        REAL_LITERAL(1.0) / std::max(cellSize, REAL_EPSILON);

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 0; j < ny; ++j) {
      for (int i = 1; i < nx; ++i) {
        const labeltype left = fields->Label(i, j + 1);
        const labeltype right = fields->Label(i + 1, j + 1);
        const bool leftFluid = IS_FLUID(left);
        const bool rightFluid = IS_FLUID(right);
        if (leftFluid == rightFluid)
          continue;
        if (IS_SOLID(left))
          continue;
        if (IS_SOLID(right))
          continue;
        const varType kappa_face =
            std::clamp(REAL_LITERAL(0.5) * (fields->kappa->Get(i - 1, j) +
                                            fields->kappa->Get(i, j)),
                       -maxAbsKappa, maxAbsKappa);
        fields->interface_u->Set(i, j, -gamma * kappa_face);
      }
    }

    OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
    for (int j = 1; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        const labeltype labelBot = fields->Label(i + 1, j);
        const labeltype labelTop = fields->Label(i + 1, j + 1);
        const bool botFluid = IS_FLUID(labelBot);
        const bool topFluid = IS_FLUID(labelTop);
        if (botFluid == topFluid)
          continue;
        if (IS_SOLID(labelBot))
          continue;
        if (IS_SOLID(labelTop))
          continue;
        const varType kappa_face =
            std::clamp(REAL_LITERAL(0.5) * (fields->kappa->Get(i, j - 1) +
                                            fields->kappa->Get(i, j)),
                       -maxAbsKappa, maxAbsKappa);
        fields->interface_v->Set(i, j, -gamma * kappa_face);
      }
    }
}

void PIC::LaplacePressure() {
  UpdatePhiFromParticles();
  ComputeSurfaceTensionOnFaces();
}
