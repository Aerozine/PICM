#include "PIC.hpp"

// Advection of all particles based on RK2.

void PIC::AdvectParticles() {

  int ppcx = particles->ppcx;
  int ppcy = particles->ppcy;

  varType xmid, ymid;
  varType umid, vmid;

  varType x1, y1;

  for (int icell = 0; icell < nx; icell++) {
    for (int jcell = 0; jcell < ny; jcell++) {

      for (int a = 0; a < ppcx; a++) {
        for (int b = 0; b < ppcy; b++) {
          int ip = icell * ppcx + a;
          int jp = jcell * ppcy + b;

          varType x0 = particles->GetX(ip, jp); // pour ne pas refaire ça 
          varType y0 = particles->GetY(ip, jp); // plusieurs fois on peut 
          varType u0 = particles->GetU(ip, jp); // fusionner ProjectGridOnParticles
          varType v0 = particles->GetV(ip, jp); // avec AdvectParticles

          xmid = x0 + 0.5 * dt * u0; 
          ymid = y0 + 0.5 * dt * v0; 

          umid = interpolateU(xmid, ymid);
          vmid = interpolateV(xmid, ymid);

          x1 = x0 + dt * umid;
          y1 = y0 + dt * vmid;

          if (x1 < 0.0 || x1 > dx * nx) {
            particles->SetDead(ip, jp, true); 
          } else if (y1 < 0.0 || y1 > dy * ny) {
            particles->SetDead(ip, jp, true); 
          } else {
            particles->SetX(ip, jp, x1);
            particles->SetY(ip, jp, y1);
          }
       }
      }
    }
  }
}

// Bilinear interpolation

varType PIC::interpolateU(const varType x, const varType y) const {
  const varType i_real = x / dx;
  const varType j_real = y / dy - REAL_LITERAL(0.5);

  int i = static_cast<int>(std::floor(i_real));
  int j = static_cast<int>(std::floor(j_real));

  const varType fx = i_real - static_cast<varType>(i);
  const varType fy = j_real - static_cast<varType>(j);

  i = std::clamp(i, 0, fields->u.nx - 2);
  j = std::clamp(j, 0, fields->u.ny - 2);

  const varType u00 = fields->u.Get(i, j);
  const varType u10 = fields->u.Get(i + 1, j);
  const varType u01 = fields->u.Get(i, j + 1);
  const varType u11 = fields->u.Get(i + 1, j + 1);

  return (REAL_LITERAL(1.0) - fy) *
             ((REAL_LITERAL(1.0) - fx) * u00 + fx * u10) +
         fy * ((REAL_LITERAL(1.0) - fx) * u01 + fx * u11);
}

varType PIC::interpolateV(const varType x, const varType y) const {
  const varType i_real = x / dx - REAL_LITERAL(0.5);
  const varType j_real = y / dy;

  int i = static_cast<int>(std::floor(i_real));
  int j = static_cast<int>(std::floor(j_real));

  const varType fx = i_real - static_cast<varType>(i);
  const varType fy = j_real - static_cast<varType>(j);

  i = std::clamp(i, 0, fields->v.nx - 2);
  j = std::clamp(j, 0, fields->v.ny - 2);

  const varType v00 = fields->v.Get(i, j);
  const varType v10 = fields->v.Get(i + 1, j);
  const varType v01 = fields->v.Get(i, j + 1);
  const varType v11 = fields->v.Get(i + 1, j + 1);

  return (REAL_LITERAL(1.0) - fy) *
             ((REAL_LITERAL(1.0) - fx) * v00 + fx * v10) +
         fy * ((REAL_LITERAL(1.0) - fx) * v01 + fx * v11);
}
