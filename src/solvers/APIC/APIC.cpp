#include "APIC.hpp"

APIC::APIC(Parameters &params) : PIC(params) {}

void APIC::ProjectGridOnParticles() {
  OMP_PRAGMA(omp parallel for)
  for (int idx = 0; idx < particles->size(); ++idx) {
    const varType xp = particles->GetX(idx);
    const varType yp = particles->GetY(idx);

    varType uNew = varType(0);
    varType vNew = varType(0);

    varType bu0 = varType(0);
    varType bu1 = varType(0);

    varType Du00 = varType(0);
    varType Du01 = varType(0);
    varType Du11 = varType(0);

    {
      const varType xg = xp / dx;
      const varType yg = yp / dy - varType(0.5);

      const int i0 = static_cast<int>(std::floor(xg));
      const int j0 = static_cast<int>(std::floor(yg));
      const int radius = params.kernelOrder;

      for (int dj = -radius; dj <= radius; ++dj) {
        for (int di = -radius; di <= radius; ++di) {
          const int i = i0 + di;
          const int j = j0 + dj;

          if (i < 0 || i >= fields->u.nx || j < 0 || j >= fields->u.ny)
            continue;

          const varType w = hat(xg - varType(i)) * hat(yg - varType(j));
          if (w <= varType(0))
            continue;

          const varType uFace = fields->u.Get(i, j);

          const varType xFace = varType(i) * dx;
          const varType yFace = (varType(j) + varType(0.5)) * dy;

          const varType ox = xFace - xp;
          const varType oy = yFace - yp;

          uNew += w * uFace;

          bu0 += w * uFace * ox;
          bu1 += w * uFace * oy;

          Du00 += w * ox * ox;
          Du01 += w * ox * oy;
          Du11 += w * oy * oy;
        }
      }
    }

    varType bv0 = varType(0);
    varType bv1 = varType(0);

    varType Dv00 = varType(0);
    varType Dv01 = varType(0);
    varType Dv11 = varType(0);

    {
      const varType xg = xp / dx - varType(0.5);
      const varType yg = yp / dy;

      const int i0 = static_cast<int>(std::floor(xg));
      const int j0 = static_cast<int>(std::floor(yg));
      const int radius = params.kernelOrder;

      for (int dj = -radius; dj <= radius; ++dj) {
        for (int di = -radius; di <= radius; ++di) {
          const int i = i0 + di;
          const int j = j0 + dj;

          if (i < 0 || i >= fields->v.nx || j < 0 || j >= fields->v.ny)
            continue;

          const varType w = hat(xg - varType(i)) * hat(yg - varType(j));
          if (w <= varType(0))
            continue;

          const varType vFace = fields->v.Get(i, j);

          const varType xFace = (varType(i) + varType(0.5)) * dx;
          const varType yFace = varType(j) * dy;

          const varType ox = xFace - xp;
          const varType oy = yFace - yp;

          vNew += w * vFace;

          bv0 += w * vFace * ox;
          bv1 += w * vFace * oy;

          Dv00 += w * ox * ox;
          Dv01 += w * ox * oy;
          Dv11 += w * oy * oy;
        }
      }
    }

    // ----- C_u = b_u D_u^{-1} -----
    varType cuX = varType(0);
    varType cuY = varType(0);

    const varType detU = Du00 * Du11 - Du01 * Du01;
    if (detU > REAL_EPSILON) {
      const varType inv00 =  Du11 / detU;
      const varType inv01 = -Du01 / detU;
      const varType inv11 =  Du00 / detU;

      cuX = bu0 * inv00 + bu1 * inv01;
      cuY = bu0 * inv01 + bu1 * inv11;
    } else {
      cuX = varType(0);
      cuY = varType(0);
    }

    // ----- C_v = b_v D_v^{-1} -----
    varType cvX = varType(0);
    varType cvY = varType(0);

    const varType detV = Dv00 * Dv11 - Dv01 * Dv01;
    if (detV > REAL_EPSILON) {
      const varType inv00 =  Dv11 / detV;
      const varType inv01 = -Dv01 / detV;
      const varType inv11 =  Dv00 / detV;

      cvX = bv0 * inv00 + bv1 * inv01;
      cvY = bv0 * inv01 + bv1 * inv11;
    } else {
      cvX = varType(0);
      cvY = varType(0);
    }

    particles->SetU(idx, uNew);
    particles->SetV(idx, vNew);
    particles->SetCu(idx, cuX, cuY);
    particles->SetCv(idx, cvX, cvY);
  }
}

void APIC::ScatterToGridAPIC(varType xg, varType yg, varType xp, varType yp,
                         varType baseVal, varType cX, varType cY, 
                         varType faceOffsetX, varType faceOffsetY,
                         Grid2D &sum, Grid2D &weight, int imax, int jmax) {
  const int i0 = static_cast<int>(std::floor(xg));
  const int j0 = static_cast<int>(std::floor(yg));

  const int radius = params.kernelOrder;

  for (int dj = -radius; dj <= radius; ++dj) {
    for (int di = -radius; di <= radius; ++di) {
      const int i = i0 + di;
      const int j = j0 + dj;

      if (i < 0 || i >= imax || j < 0 || j >= jmax)
        continue;

      const varType w = hat(xg - varType(i)) * hat(yg - varType(j));

      if (w <= varType(0))
        continue;

      const varType xFace = (varType(i) + faceOffsetX) * dx;
      const varType yFace = (varType(j) + faceOffsetY) * dy;

      const varType ox = xFace - xp;
      const varType oy = yFace - yp;

      const varType affineVal = baseVal + cX * ox + cY * oy;

      sum.Set(i, j, sum.Get(i, j) + w * affineVal);
      weight.Set(i, j, weight.Get(i, j) + w);
    }
  }
}

void APIC::ProjectParticleOnMAC(int idx) {
  varType x = particles->GetX(idx);
  varType y = particles->GetY(idx);
  varType up = particles->GetU(idx);
  varType vp = particles->GetV(idx);
  varType cuX = particles->GetCuX(idx);
  varType cuY = particles->GetCuY(idx);
  varType cvX = particles->GetCvX(idx);
  varType cvY = particles->GetCvY(idx);

  ScatterToGridAPIC(x / dx, y / dy - varType(0.5), x, y, up, cuX, cuY,
                  varType(0), varType(0.5), *fields->u_sum,
                  *fields->u_weight, fields->u.nx, fields->u.ny);

  ScatterToGridAPIC(x / dx - varType(0.5), y / dy, x, y, vp, cvX, cvY,
                  varType(0.5), varType(0), *fields->v_sum,
                  *fields->v_weight, fields->v.nx, fields->v.ny);
}
