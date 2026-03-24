#include "PIC.hpp"
#include <iostream>

varType PIC::GetW() {
  return static_cast<varType>(particles->ppcx * particles->ppcy);
}

varType PIC::hat(varType r) {
  if (r >= varType(0) && r <= varType(1))
    return varType(1) - r;
  else if (r >= varType(-1) && r < varType(0))
    return varType(1) + r;
  else
    return varType(0);
}

void PIC::ProjectOneParticleOnMAC(varType x, varType y, varType up,
                                  varType vp) {
  // u faces: staggered at (i·dx, (j+0.5)·dy)
  varType xu = x / dx;
  varType yu = y / dy - varType(0.5);
  int i0 = static_cast<int>(std::floor(xu));
  int j0 = static_cast<int>(std::floor(yu));

  for (int dj = 0; dj <= 1; ++dj) {
    for (int di = 0; di <= 1; ++di) {
      int i = i0 + di, j = j0 + dj;
      if (i < 0 || i >= nx + 1 || j < 0 || j >= ny)
        continue;
      varType k = hat(xu - varType(i)) * hat(yu - varType(j));
      if (k > varType(0)) {
        fields->u_sum.Set(i, j, fields->u_sum.Get(i, j) + k * up);
        fields->u_weight.Set(i, j, fields->u_weight.Get(i, j) + k);
      }
    }
  }

  // v faces: staggered at ((i+0.5)·dx, j·dy)
  varType xv = x / dx - varType(0.5);
  varType yv = y / dy;
  i0 = static_cast<int>(std::floor(xv));
  j0 = static_cast<int>(std::floor(yv));

  for (int dj = 0; dj <= 1; ++dj) {
    for (int di = 0; di <= 1; ++di) {
      int i = i0 + di, j = j0 + dj;
      if (i < 0 || i >= nx || j < 0 || j >= ny + 1)
        continue;
      varType k = hat(xv - varType(i)) * hat(yv - varType(j));
      if (k > varType(0)) {
        fields->v_sum.Set(i, j, fields->v_sum.Get(i, j) + k * vp);
        fields->v_weight.Set(i, j, fields->v_weight.Get(i, j) + k);
      }
    }
  }
}

void PIC::ProjectParticlesOnGrid(std::string kernel) {
  if (kernel != "hat") {
    std::cerr << "[PIC] Unknown P2G kernel '" << kernel << "'.\n";
    return;
  }

  // Zero accumulation buffers.
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx + 1; i++) {
      fields->u_sum.Set(i, j, varType(0));
      fields->u_weight.Set(i, j, varType(0));
    }
  for (int j = 0; j < ny + 1; j++)
    for (int i = 0; i < nx; i++) {
      fields->v_sum.Set(i, j, varType(0));
      fields->v_weight.Set(i, j, varType(0));
    }

  // Accumulate — iterate flat over all allocated slots.
  const int cap = particles->capacity;
  for (int idx = 0; idx < cap; ++idx) {
    if (particles->IsDead(idx))
      continue;
    ProjectOneParticleOnMAC(particles->GetX(idx), particles->GetY(idx),
                            particles->GetU(idx), particles->GetV(idx));
  }

  // Normalize and write back, preserving BC and solid faces.
  // u faces: face (i,j) is between cells (i-1,j) and (i,j).
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx + 1; i++) {
      bool isSolid = false, isBC = false;
      if (i > 0) {
        auto l = fields->Label(i - 1, j);
        if (l & Fields2D::SOLID)
          isSolid = true;
        if (l & Fields2D::BC_U)
          isBC = true;
      }
      if (i < nx) {
        auto l = fields->Label(i, j);
        if (l & Fields2D::SOLID)
          isSolid = true;
        if (l & Fields2D::BC_U)
          isBC = true;
      }
      if (isSolid)
        fields->u.Set(i, j, fields->usolid);
      else if (isBC)
        ; // keep value written by applyToFields
      else if (fields->u_weight.Get(i, j) > varType(1e-12))
        fields->u.Set(i, j,
                      fields->u_sum.Get(i, j) / fields->u_weight.Get(i, j));
      else
        fields->u.Set(i, j, varType(0));
    }
  }

  // v faces: face (i,j) is between cells (i,j-1) and (i,j).
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx; i++) {
      bool isSolid = false, isBC = false;
      if (j > 0) {
        auto l = fields->Label(i, j - 1);
        if (l & Fields2D::SOLID)
          isSolid = true;
        if (l & Fields2D::BC_V)
          isBC = true;
      }
      if (j < ny) {
        auto l = fields->Label(i, j);
        if (l & Fields2D::SOLID)
          isSolid = true;
        if (l & Fields2D::BC_V)
          isBC = true;
      }
      if (isSolid)
        fields->v.Set(i, j, fields->usolid);
      else if (isBC)
        ;
      else if (fields->v_weight.Get(i, j) > varType(1e-12))
        fields->v.Set(i, j,
                      fields->v_sum.Get(i, j) / fields->v_weight.Get(i, j));
      else
        fields->v.Set(i, j, varType(0));
    }
  }
}

void PIC::ProjectGridOnParticles() {
  const int cap = particles->capacity;
  for (int idx = 0; idx < cap; ++idx) {
    if (particles->IsDead(idx))
      continue;
    particles->SetU(idx,
                    interpolateU(particles->GetX(idx), particles->GetY(idx)));
    particles->SetV(idx,
                    interpolateV(particles->GetX(idx), particles->GetY(idx)));
  }
}
