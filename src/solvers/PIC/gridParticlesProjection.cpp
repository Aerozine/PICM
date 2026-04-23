#include "PIC.hpp"
#include <iostream>


void PIC::ScatterToGrid(varType xg, varType yg, varType val, Grid2D &sum,
                        Grid2D &weight, int imax, int jmax) {
  int i0 = static_cast<int>(std::floor(xg));
  int j0 = static_cast<int>(std::floor(yg));

#ifdef USE_SPEED
  constexpr int8_t radius=1;
#else
  constexpr int8_t radius=2;
#endif

  //int radius = params.kernelOrder;
  OMP_PRAGMA(omp parallel for collapse(2))
  for (int8_t dj = -radius; dj <= radius; ++dj) {
    for (int di = -radius; di <= radius; ++di) {
      int i = i0 + di, j = j0 + dj;
      if (i < 0 || i >= imax || j < 0 || j >= jmax)
        continue;
      varType k = hat(xg - varType(i)) * hat(yg - varType(j));
      if (k > varType(0)) {
    OMP_PRAGMA(omp critical)
        {
          sum.Set(i, j, sum.Get(i, j) + k * val);
          weight.Set(i, j, weight.Get(i, j) + k);
        }
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
  fields->u_sum->reset();
  fields->u_weight->reset();

  fields->v_sum->reset();
  fields->v_weight->reset();

  OMP_PRAGMA(omp parallel for collapse(2))
  for (int idx = 0; idx < particles->size(); ++idx) {
    ProjectParticleOnMAC(idx);
  }

  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < fields->u.ny; j++) {
    for (int i = 0; i < fields->u.nx; i++) {
      labeltype left = fields->Label(i, j + 1);
      labeltype right = fields->Label(i + 1, j + 1);

      if (IS_SOLID(left) || IS_SOLID(right)) {
        fields->u.Set(i, j, varType(0));
        continue;
      } if (IS_BC_U(left)) {
        continue;
      }

      if (fields->u_weight->Get(i, j) >= REAL_EPSILON) {
        varType uNew = fields->u_sum->Get(i, j) / fields->u_weight->Get(i, j);
        fields->u.Set(i, j, uNew);
      } else if (!IS_BC_U(left))
          fields->u.Set(i, j, varType(0));
    }
  }

  OMP_PRAGMA(omp parallel for collapse(2))
  for (int j = 0; j < fields->v.ny; j++) {
    for (int i = 0; i < fields->v.nx; i++) {
      labeltype bottom = fields->Label(i + 1, j);
      labeltype top = fields->Label(i + 1, j + 1);

      if (IS_SOLID(bottom) || IS_SOLID(top)) {
        fields->v.Set(i, j, varType(0));
        continue;
      } if (IS_BC_V(bottom)) {
        continue;
      }
      if (fields->v_weight->Get(i, j) >= REAL_EPSILON) {
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
  varType x = particles->GetX(idx);
  varType y = particles->GetY(idx);
  particles->SetU(idx,
    fields->u.interpolate(x,y,dx,dy,0));
  particles->SetV(idx,
    fields->v.interpolate(x,y,dx,dy,1) - dt*params.gravity);
}
}

// same function as above in PIC but needed for FLIP
// because ProjectGridOnParticles is virtual (rewritten for FLIP)
/*void PIC::ProjectBCOnParticles() {
OMP_PRAGMA(omp parallel for)
for (int idx = 0; idx < particles->size(); ++idx) {
  particles->SetU(
      idx,
    u.interpolate(x, y, dx, dy, 0);
      interpolateU(fields->u, particles->GetX(idx), particles->GetY(idx)));
  particles->SetV(
      idx, interpolateV(fields->v, particles->GetX(idx), particles->GetY(idx)));
}
}

*/