#include "PIC.hpp"

void PIC::AdvectParticles() {
  const varType xMax = dx * nx;
  const varType yMax = dy * ny;
  const int cap = particles->capacity;

  for (int idx = 0; idx < cap; ++idx) {
    if (particles->IsDead(idx))
      continue;

    varType x0 = particles->GetX(idx);
    varType y0 = particles->GetY(idx);
    varType u0 = particles->GetU(idx);
    varType v0 = particles->GetV(idx);

    // RK2 mid-point
    varType xmid = x0 + varType(0.5) * dt * u0;
    varType ymid = y0 + varType(0.5) * dt * v0;
    varType umid = interpolateU(xmid, ymid);
    varType vmid = interpolateV(xmid, ymid);

    varType x1 = x0 + dt * umid;
    varType y1 = y0 + dt * vmid;

    // Kill particles that left the domain.
    if (x1 < varType(0) || x1 >= xMax || y1 < varType(0) || y1 >= yMax) {
      particles->SetDead(idx, true);
      deadSlots->push(idx);
      continue;
    }

    int i1 = std::clamp(static_cast<int>(std::floor(x1 / dx)), 0, nx - 1);
    int j1 = std::clamp(static_cast<int>(std::floor(y1 / dy)), 0, ny - 1);

    if (fields->Label(i1, j1) & Fields2D::SOLID) {
      particles->SetDead(idx, true);
      deadSlots->push(idx);
      continue;
    } else {
      particles->SetX(idx, x1);
      particles->SetY(idx, y1);
    }
  }
}

varType PIC::interpolateU(const varType x, const varType y) const {
  const varType i_real = x / dx;
  const varType j_real = y / dy - REAL_LITERAL(0.5);
  int i = static_cast<int>(std::floor(i_real));
  int j = static_cast<int>(std::floor(j_real));
  const varType fx = i_real - varType(i);
  const varType fy = j_real - varType(j);
  i = std::clamp(i, 0, fields->u.nx - 2);
  j = std::clamp(j, 0, fields->u.ny - 2);
  return (1 - fy) *
             ((1 - fx) * fields->u.Get(i, j) + fx * fields->u.Get(i + 1, j)) +
         fy * ((1 - fx) * fields->u.Get(i, j + 1) +
               fx * fields->u.Get(i + 1, j + 1));
}

varType PIC::interpolateV(const varType x, const varType y) const {
  const varType i_real = x / dx - REAL_LITERAL(0.5);
  const varType j_real = y / dy;
  int i = static_cast<int>(std::floor(i_real));
  int j = static_cast<int>(std::floor(j_real));
  const varType fx = i_real - varType(i);
  const varType fy = j_real - varType(j);
  i = std::clamp(i, 0, fields->v.nx - 2);
  j = std::clamp(j, 0, fields->v.ny - 2);
  return (1 - fy) *
             ((1 - fx) * fields->v.Get(i, j) + fx * fields->v.Get(i + 1, j)) +
         fy * ((1 - fx) * fields->v.Get(i, j + 1) +
               fx * fields->v.Get(i + 1, j + 1));
}
