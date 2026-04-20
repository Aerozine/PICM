#include "APIC.hpp"

APIC::APIC(Parameters &params) : PIC(params) {}

// todo needs to be fixed for each type of hat
varType APIC::dhat(varType r) {
  if (varType(-1.5) <= r && r < varType(-0.5))
    return r + static_cast<varType>(1.5);
  if (varType(-0.5) <= r && r < varType(0.5))
    return -static_cast<varType>(2.0) * r;
  if (varType(0.5) <= r && r < varType(1.5))
    return r - static_cast<varType>(1.5);
  return static_cast<varType>(0);
}

Vec2 APIC::gradWeightU(const int i, const int j, const varType xp,
                       const varType yp) const {
  const varType rx = xp / dx - static_cast<varType>(i);
  const varType ry =
      yp / dy - static_cast<varType>(0.5) - static_cast<varType>(j);

  Vec2 g;
  g.x = (dhat(rx) * hat(ry)) / dx;
  g.y = (hat(rx) * dhat(ry)) / dy;
  return g;
}

Vec2 APIC::gradWeightV(int i, int j, varType xp, varType yp) const {
  const varType rx =
      xp / dx - static_cast<varType>(0.5) - static_cast<varType>(i);
  const varType ry = yp / dy - static_cast<varType>(j);

  Vec2 g;
  g.x = dhat(rx) * hat(ry) / dx;
  g.y = hat(rx) * dhat(ry) / dy;
  return g;
}

void APIC::ProjectGridOnParticles() {
    OMP_PRAGMA(omp parallel for)
    for (int idx = 0; idx < particles->size(); ++idx) {
      const varType xp = particles->GetX(idx);
      const varType yp = particles->GetY(idx);

      // 1) vitesse particule : comme PIC
      const varType uNew = interpolateU(fields->u, xp, yp);
      const varType vNew = interpolateV(fields->v, xp, yp);

      // 2) reconstruction de cu
      Vec2 cu{static_cast<varType>(0), static_cast<varType>(0)};
      {
        const varType i_real = xp / dx;
        const varType j_real = yp / dy - static_cast<varType>(0.5);
        const int i0 = static_cast<int>(std::floor(i_real));
        const int j0 = static_cast<int>(std::floor(j_real));

        for (int dj = -1; dj <= 1; ++dj) {
          for (int di = -1; di <= 1; ++di) {
            const int i = i0 + di;
            const int j = j0 + dj;
            if (i < 0 || i >= fields->u.nx || j < 0 || j >= fields->u.ny)
              continue;

            Vec2 gw = gradWeightU(i, j, xp, yp);
            const varType uFace = fields->u.Get(i, j);

            cu.x += gw.x * uFace;
            cu.y += gw.y * uFace;
          }
        }
      }

      // 3) reconstruction de cv
      Vec2 cv{varType(0), varType(0)};
      {
        const varType i_real = xp / dx - varType(0.5);
        const varType j_real = yp / dy;
        const int i0 = static_cast<int>(std::floor(i_real));
        const int j0 = static_cast<int>(std::floor(j_real));

        for (int dj = -1; dj <= 1; ++dj) {
          for (int di = -1; di <= 1; ++di) {
            const int i = i0 + di;
            const int j = j0 + dj;
            if (i < 0 || i >= fields->v.nx || j < 0 || j >= fields->v.ny)
              continue;

            Vec2 gw = gradWeightV(i, j, xp, yp);
            const varType vFace = fields->v.Get(i, j);

            cv.x += gw.x * vFace;
            cv.y += gw.y * vFace;
          }
        }
      }

      particles->SetU(idx, uNew);
      particles->SetV(idx, vNew);
      particles->SetCu(idx, cu.x, cu.y);
      particles->SetCv(idx, cv.x, cv.y);
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

  // u faces
  {
    const varType xg = x / dx;
    const varType yg = y / dy - varType(0.5);

    const int i0 = static_cast<int>(std::floor(xg));
    const int j0 = static_cast<int>(std::floor(yg));

    for (int dj = -1; dj <= 1; ++dj) {
      for (int di = -1; di <= 1; ++di) {
        const int i = i0 + di;
        const int j = j0 + dj;
        if (i < 0 || i >= fields->u.nx || j < 0 || j >= fields->u.ny)
          continue;

        const varType w = hat(xg - varType(i)) * hat(yg - varType(j));
        if (w <= varType(0))
          continue;

        const varType xFace = varType(i) * dx;
        const varType yFace = (varType(j) + varType(0.5)) * dy;

        const varType dux = xFace - x;
        const varType duy = yFace - y;

        const varType uAff = up + cuX * dux + cuY * duy;

        fields->u_sum->Set(i, j, fields->u_sum->Get(i, j) + w * uAff);
        fields->u_weight->Set(i, j, fields->u_weight->Get(i, j) + w);
      }
    }
  }

  // v faces
  {
    const varType xg = x / dx - varType(0.5);
    const varType yg = y / dy;

    const int i0 = static_cast<int>(std::floor(xg));
    const int j0 = static_cast<int>(std::floor(yg));

    for (int dj = -1; dj <= 1; ++dj) {
      for (int di = -1; di <= 1; ++di) {
        const int i = i0 + di;
        const int j = j0 + dj;
        if (i < 0 || i >= fields->v.nx || j < 0 || j >= fields->v.ny)
          continue;

        const varType w = hat(xg - varType(i)) * hat(yg - varType(j));
        if (w <= varType(0))
          continue;

        const varType xFace = (varType(i) + varType(0.5)) * dx;
        const varType yFace = varType(j) * dy;

        const varType dvx = xFace - x;
        const varType dvy = yFace - y;

        const varType vAff = vp + cvX * dvx + cvY * dvy;

        fields->v_sum->Set(i, j, fields->v_sum->Get(i, j) + w * vAff);
        fields->v_weight->Set(i, j, fields->v_weight->Get(i, j) + w);
      }
    }
  }
}