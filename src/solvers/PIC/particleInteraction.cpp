#include "PIC.hpp"

void PIC::particleInteraction() {
  const varType cellSize = std::min(dx, dy);
  const varType nominalParticleSpacing =
      std::min(dx / static_cast<varType>(std::max(1, params.ppcx)),
               dy / static_cast<varType>(std::max(1, params.ppcy)));
  const varType particleRadius =
      params.particleRadius > REAL_EPSILON
          ? params.particleRadius
          : REAL_LITERAL(0.5) * nominalParticleSpacing;
  const varType interactionRadius =
      std::max(REAL_LITERAL(6.0) * particleRadius,
               REAL_LITERAL(2.0) * cellSize);
  const int neighborCellRadius =
      std::max(1, static_cast<int>(std::ceil(interactionRadius / cellSize)));
  const int interfaceBandCells = 2;
  const varType kernelExponent =
      std::max(params.interactionExponent, REAL_LITERAL(1.0));
  const varType minDt = std::max(dt, REAL_LITERAL(1e-12));
  const varType maxVelocityKick = REAL_LITERAL(0.25) * cellSize / minDt;

  OMP_PRAGMA(omp parallel for collapse(2) schedule(dynamic,1))
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const labeltype center = fields->Label(i + 1, j + 1);
      if (!IS_FLUID(center))
        continue;

      int ci_lo = std::max(0, i - interfaceBandCells);
      int ci_hi = std::min(nx - 1, i + interfaceBandCells);
      int cj_lo = std::max(0, j - interfaceBandCells);
      int cj_hi = std::min(ny - 1, j + interfaceBandCells);

      bool nearAir = false;
      for (int cj = cj_lo; cj <= cj_hi && !nearAir; ++cj) {
        for (int ci = ci_lo; ci <= ci_hi; ++ci) {
          const labeltype left = fields->Label(ci, cj + 1);
          const labeltype right = fields->Label(ci + 2, cj + 1);
          const labeltype up = fields->Label(ci + 1, cj);
          const labeltype down = fields->Label(ci + 1, cj + 2);
          nearAir = IS_AIR(left) || IS_AIR(right) || IS_AIR(up) || IS_AIR(down);
          if (nearAir)
            break;
        }
      }
      if (!nearAir)
        continue;

      Particles &centerCell = (*cloud)(i, j);
      const int particleCount = centerCell.size();
      if (particleCount == 0)
        continue;

      ci_lo = std::max(0, i - neighborCellRadius);
      ci_hi = std::min(nx - 1, i + neighborCellRadius);
      cj_lo = std::max(0, j - neighborCellRadius);
      cj_hi = std::min(ny - 1, j + neighborCellRadius);

      for (int p = 0; p < particleCount; ++p) {
        varType weightedDx = 0.0;
        varType weightedDy = 0.0;
        varType weightSum = 0.0;

        const varType xp = centerCell.GetX(p);
        const varType yp = centerCell.GetY(p);

        for (int cj = cj_lo; cj <= cj_hi; ++cj) {
          for (int ci = ci_lo; ci <= ci_hi; ++ci) {
            if (!IS_FLUID(fields->Label(ci + 1, cj + 1)))
              continue;

            const Particles &neighborCell = (*cloud)(ci, cj);
            const int neighborParticleCount = neighborCell.size();

            for (int q = 0; q < neighborParticleCount; ++q) {
              if (ci == i && cj == j && q == p)
                continue;

              const varType dxpq = neighborCell.GetX(q) - xp;
              const varType dypq = neighborCell.GetY(q) - yp;
              const varType r2 = dxpq * dxpq + dypq * dypq;
              if (r2 < REAL_EPSILON)
                continue;

              const varType r = std::sqrt(r2);
              if (r >= interactionRadius)
                continue;

              const varType qn = 1.0 - (r / interactionRadius);
              const varType weight = std::pow(qn, kernelExponent);

              weightedDx += weight * dxpq;
              weightedDy += weight * dypq;
              weightSum += weight;
            }
          }
        }

        if (weightSum <= REAL_EPSILON)
          continue;

        // The centroid mismatch is close to zero in the bulk and points
        // inward near a free surface, which avoids singular 1/r forces.
        varType du = dt * params.interactionStiffness * weightedDx /
                     (weightSum * interactionRadius);
        varType dv = dt * params.interactionStiffness * weightedDy /
                     (weightSum * interactionRadius);

        const varType kickNorm = std::hypot(du, dv);
        if (kickNorm > maxVelocityKick) {
          const varType scale = maxVelocityKick / kickNorm;
          du *= scale;
          dv *= scale;
        }

        centerCell.SetU(p, centerCell.GetU(p) + du);
        centerCell.SetV(p, centerCell.GetV(p) + dv);
      }
    }
  }
}
