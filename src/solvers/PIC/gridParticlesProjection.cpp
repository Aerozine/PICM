#include "PIC.hpp"
#include <iostream>

varType PIC::GetW() {
  return static_cast<varType>(particles->ppcx * particles->ppcy);
}
// @todo handle different hat correctly
varType hat1(varType r) {
  if (r >= varType(0) && r <= varType(1))
    return varType(1) - r;
  else if (r >= varType(-1) && r < varType(0))
    return varType(1) + r;
  else
    return varType(0);
}
// h2
varType PIC::hat(varType r) {
  if (varType(-1.5) <= r && r < varType(-0.5))
    return varType(0.5) * (r + varType(3.0 / 2.0)) * (r + varType(3.0 / 2.0));
  if (varType(-0.5) <= r && r < varType(0.5))
    return varType(0.75) - r * r;
  if (varType(0.5) <= r && r < varType(1.5))
    return varType(0.5) * (varType(3.0 / 2.0) - r) * (varType(3.0 / 2.0) - r);
  return varType(0);
}
void PIC::ScatterToGrid(varType xg, varType yg, varType val, Grid2D &sum,
                        Grid2D &weight, int imax, int jmax) {
  int i0 = static_cast<int>(std::floor(xg));
  int j0 = static_cast<int>(std::floor(yg));
  // @todo dynamic to hat size ! for h1 2 , h2 3
  for (int dj = -1; dj <= 1; ++dj) {
    for (int di = -1; di <= 1; ++di) {
      int i = i0 + di, j = j0 + dj;
      if (i < 0 || i >= imax || j < 0 || j >= jmax)
        continue;
      varType k = hat(xg - varType(i)) * hat(yg - varType(j));
      if (k > varType(0)) {
        sum.Set(i, j, sum.Get(i, j) + k * val);
        weight.Set(i, j, weight.Get(i, j) + k);
      }
    }
  }
}
void PIC::ProjectParticleOnMAC(varType x, varType y, varType up, varType vp) {
  ScatterToGrid(x / dx, y / dy - varType(0.5), up, fields->u_sum,
                fields->u_weight, nx + 1, ny);
  ScatterToGrid(x / dx - varType(0.5), y / dy, vp, fields->v_sum,
                fields->v_weight, nx, ny + 1);
}

void PIC::ProjectParticlesOnGrid(std::string kernel) {
  if (kernel != "hat") {
    std::cerr << "[PIC] Unknown P2G kernel '" << kernel << "'.\n";
    return;
  }

  // Weights to zero.
OMP_PRAGMA(omp parallel for collapse(2))
for (int j = 0; j < fields->u.ny; j++)
  for (int i = 0; i < fields->u.nx; i++) {
    fields->u_sum.Set(i, j, varType(0));
    fields->u_weight.Set(i, j, varType(0));
  }

OMP_PRAGMA(omp parallel for collapse(2))
for (int j = 0; j < fields->v.ny; j++)
  for (int i = 0; i < fields->v.nx; i++) {
    fields->v_sum.Set(i, j, varType(0));
    fields->v_weight.Set(i, j, varType(0));
  }
// will messed up the code due to race condition and other funny things
// @todo must be implemented directly inside the project particle on MAC
// OMP_PRAGMA(omp parallel for)
for (int idx = 0; idx < particles->size(); ++idx) {
  ProjectParticleOnMAC(particles->GetX(idx), particles->GetY(idx),
                       particles->GetU(idx), particles->GetV(idx));
}

// Normalize u faces.
// u(i,j) sits between Label(i, j+1) [left] and Label(i+1, j+1) [right].
// This matches SemiLagrangian::Advect() u-loop convention.
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < fields->u.ny; j++) {
    for (int i = 0; i < fields->u.nx; i++) {
      labeltype left = fields->Label(i, j + 1);
      labeltype right = fields->Label(i + 1, j + 1);

      if (IS_SOLID(left) || IS_SOLID(right)) {
        fields->u.Set(i, j, varType(0));
        continue;
      }
      if (IS_BC_U(left)) {
        // Boundary value already set — do not overwrite
        continue;
      }
      //@todo should be illegal to do 1e-12
      if (fields->u_weight.Get(i, j) > varType(1e-12))
        fields->u.Set(i, j,
                      fields->u_sum.Get(i, j) / fields->u_weight.Get(i, j));
      else if (!IS_BC_U(left))
        fields->u.Set(i, j, varType(0));
    }
  }

  // Normalize v faces.
  // v(i,j) sits between Label(i+1, j) [bottom] and Label(i+1, j+1) [top].
  // This matches SemiLagrangian::Advect() v-loop convention.
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < fields->v.ny; j++) {
    for (int i = 0; i < fields->v.nx; i++) {
      labeltype bottom = fields->Label(i + 1, j);
      labeltype top = fields->Label(i + 1, j + 1);

      if (IS_SOLID(bottom) || IS_SOLID(top)) {
        fields->v.Set(i, j, varType(0));
        continue;
      }
      if (IS_BC_V(bottom)) {
        // Boundary value already set — do not overwrite
        continue;
      }
      if (fields->v_weight.Get(i, j) > varType(1e-12)){
        varType newV = fields->v_sum.Get(i, j) / fields->v_weight.Get(i, j);
        varType oldV = fields->v.Get(i, j);
        //assert(oldV - 9.81*fields->dt - 0.0001 < newV && newV < oldV - 9.81*fields->dt + 0.0001);
        if(!(oldV - 9.81*fields->dt - 0.001 < newV && newV < oldV - 9.81*fields->dt + 0.001))
          printf("new :%f \t old: %f \t  exact :%f \n",newV, oldV,- 9.81*fields->dt );
        fields->v.Set(i, j,newV);
      }
      else if (!IS_BC_V(bottom))
        fields->v.Set(i, j, varType(0));
    }
  }
}
void PIC::ProjectGridOnParticles() {
OMP_PRAGMA(omp parallel for)
for (int idx = 0; idx < particles->size(); ++idx) {
  particles->SetU(idx,
                  interpolateU(particles->GetX(idx), particles->GetY(idx)));
  particles->SetV(idx,
                  interpolateV(particles->GetX(idx), particles->GetY(idx)));
}
}
