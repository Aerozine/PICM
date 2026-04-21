#include "PIC.hpp"
#include <iostream>

varType PIC::GetW() {
  return static_cast<varType>(particles->ppcx * particles->ppcy);
}

void PIC::ScatterToGrid(varType xg, varType yg, varType val, Grid2D &sum,
                        Grid2D &weight, int imax, int jmax) {
  int i0 = static_cast<int>(std::floor(xg));
  int j0 = static_cast<int>(std::floor(yg));

  int radius = params.kernelOrder;

  for (int dj = -radius; dj <= radius; ++dj) {
    for (int di = -radius; di <= radius; ++di) {
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
void PIC::ProjectParticleOnMAC(int idx) {
  varType x = particles->GetX(idx);
  varType y = particles->GetY(idx);
  varType up = particles->GetU(idx);
  varType vp = particles->GetV(idx);

  ScatterToGrid(x / dx, y / dy - varType(0.5), up, *fields->u_sum,
                *fields->u_weight, nx + 1, ny);
  ScatterToGrid(x / dx - varType(0.5), y / dy, vp, *fields->v_sum,
                *fields->v_weight, nx, ny + 1);
}

void PIC::ProjectParticlesOnGrid() {

  // Weights to zero.
OMP_PRAGMA(omp parallel for collapse(2))
for (int j = 0; j < fields->u.ny; j++)
  for (int i = 0; i < fields->u.nx; i++) {
    fields->u_sum->Set(i, j, varType(0));
    fields->u_weight->Set(i, j, varType(0));
  }

OMP_PRAGMA(omp parallel for collapse(2))
for (int j = 0; j < fields->v.ny; j++)
  for (int i = 0; i < fields->v.nx; i++) {
    fields->v_sum->Set(i, j, varType(0));
    fields->v_weight->Set(i, j, varType(0));
  }
// will messed up the code due to race condition and other funny things
// @todo must be implemented directly inside the project particle on MAC
// OMP_PRAGMA(omp parallel for)
for (int idx = 0; idx < particles->size(); ++idx) {
  ProjectParticleOnMAC(idx);
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
      if (fields->u_weight->Get(i, j) <=
          static_cast<varType>(100) * std::numeric_limits<varType>::epsilon()) {
        if (!IS_BC_U(left)) {
          fields->u.Set(i, j, static_cast<varType>(0));
        }
      } else {
        fields->u.Set(i, j,
                      fields->u_sum->Get(i, j) / fields->u_weight->Get(i, j));
      }
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
      if (fields->v_weight->Get(i, j) > varType(1e-12)) {
        varType newV = fields->v_sum->Get(i, j) / fields->v_weight->Get(i, j);
        fields->v.Set(i, j, newV);
      } else if (!IS_BC_V(bottom))
        fields->v.Set(i, j, varType(0));
    }
  }
}
void PIC::ProjectGridOnParticles() {
OMP_PRAGMA(omp parallel for)
for (int idx = 0; idx < particles->size(); ++idx) {
  particles->SetU(
      idx, interpolateU(fields->u, particles->GetX(idx), particles->GetY(idx)));
  particles->SetV(
      idx, interpolateV(fields->v, particles->GetX(idx), particles->GetY(idx)));
}
}

// same function as above in PIC but needed for FLIP
// because ProjectGridOnParticles is virtual (rewritten for FLIP)
void PIC::ProjectBCOnParticles() {
OMP_PRAGMA(omp parallel for)
for (int idx = 0; idx < particles->size(); ++idx) {
  particles->SetU(
      idx, interpolateU(fields->u, particles->GetX(idx), particles->GetY(idx)));
  particles->SetV(
      idx, interpolateV(fields->v, particles->GetX(idx), particles->GetY(idx)));
}
}
