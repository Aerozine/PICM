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

void PIC::ProjectParticleOnMAC(varType x, varType y, varType up, varType vp) {
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

  // Weights to zero.
  for (int j = 0; j < fields->u.ny; j++)
    for (int i = 0; i < fields->u.nx; i++) {
      fields->u_sum.Set(i, j, varType(0));
      fields->u_weight.Set(i, j, varType(0));
    }
  for (int j = 0; j < fields->v.ny; j++)
    for (int i = 0; i < fields->v.nx; i++) {
      fields->v_sum.Set(i, j, varType(0));
      fields->v_weight.Set(i, j, varType(0));
    }

  for (int idx = 0; idx < particles->size(); ++idx) {
    ProjectParticleOnMAC(particles->GetX(idx), particles->GetY(idx),
                            particles->GetU(idx), particles->GetV(idx));
  }

  // Normalize preserving BC and solid faces.
  // u faces: face (i,j) is between cells (i,j+1) and (i+1,j+1).
  for (int j = 0; j < fields->u.ny; j++) {
    for (int i = 0; i < fields->u.nx; i++) {
      bool isSolid = false, isBC = false; //, isAir = false;
      
      auto l = fields->Label(i, j + 1);
      if (l & Fields2D::SOLID)
        isSolid = true;
      if (l & Fields2D::BC_U)
        isBC = true;
      // if (l & Fields2D::AIR)
      //   isAir = true;

      auto r = fields->Label(i + 1, j + 1);
      if (r & Fields2D::SOLID)
        isSolid = true;
      if (r & Fields2D::BC_U)
        isBC = true;
      // if (r & Fields2D::AIR)
      //   isAir = true;

      if (isSolid)
        fields->u.Set(i, j, FIELD_USOLID);
      else if (isBC)
        ;
      // else if (isAir)
      //   fields->u.Set(i, j, varType(0));
      else if (fields->u_weight.Get(i, j) > varType(1e-12))
        fields->u.Set(i, j,
                      fields->u_sum.Get(i, j) / fields->u_weight.Get(i, j));
      else
        continue;
        // fields->u.Set(i, j, varType(0));
    }
  }

  // v faces: face (i,j) is between cells (i+1,j) and (i+1,j+1).
  for (int j = 0; j < fields->v.ny; j++) {
    for (int i = 0; i < fields->v.nx; i++) {
      bool isSolid = false, isBC = false; //, isAir = false;
      
      auto b = fields->Label(i + 1, j);
      if (b & Fields2D::SOLID)
        isSolid = true;
      if (b & Fields2D::BC_V)
        isBC = true;
      // if (b & Fields2D::AIR)
      //   isAir = true;
      
      auto t = fields->Label(i + 1, j + 1);
      if (t & Fields2D::SOLID)
        isSolid = true;
      if (t & Fields2D::BC_V)
        isBC = true;
      // if (t & Fields2D::AIR)
      //   isAir = true;
      
      if (isSolid)
        fields->v.Set(i, j,FIELD_USOLID);
      else if (isBC)
        ;
      // else if (isAir)
      //   fields->v.Set(i, j, varType(0));
      else if (fields->v_weight.Get(i, j) > varType(1e-12))
        fields->v.Set(i, j,
                      fields->v_sum.Get(i, j) / fields->v_weight.Get(i, j));
      else
        continue;
        // fields->v.Set(i, j, varType(0));
    }
  }
}

void PIC::ProjectGridOnParticles() {
  for (int idx = 0; idx < particles->size(); ++idx) {
    particles->SetU(idx, interpolateU(particles->GetX(idx), particles->GetY(idx)));
    particles->SetV(idx, interpolateV(particles->GetX(idx), particles->GetY(idx)));
  }
}
