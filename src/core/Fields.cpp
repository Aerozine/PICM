#include "Fields.hpp"
#include <cmath>
#include <iostream>

void Fields2D::Div() {

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      varType dudx = (u.Get(i + 1, j) - u.Get(i, j)) / dx;
      varType dvdy = (v.Get(i, j + 1) - v.Get(i, j)) / dy;

      div.Set(i, j, dudx + dvdy);
    }
  }
  return;
}

void Fields2D::VelocityNormCenterGrid() {

  for(int i = 0; i < nx; i++) {
    for(int j = 0; j < ny; j++) {

      varType x = (i * this -> dx) + dx / 2;
      varType y = (j * this -> dy) + dy / 2;

      varType uCenter = u.Interpolate(x, y, this -> dx, this -> dy, 0);
      varType vCenter = v.Interpolate(x, y, this -> dx, this -> dy, 1);
      
      varType Norm = sqrt(uCenter*uCenter + vCenter*vCenter);
      normVelocity.Set(i, j, Norm);
    }
  }
}

// FANCY CONFIG 
// simply seting the value of u as in the center
void Fields2D::InitTaylorGreen(const varType amplitude) {
    const varType Lx = static_cast<varType>(nx - 1) * dx;
    const varType Ly = static_cast<varType>(ny - 1) * dy;

    #ifdef NDEBUG
        std::cout << "Initializing Taylor-Green vortex with amplitude "
               << amplitude << std::endl;
        std::cout << "Domain: Lx = " << Lx 
                  << ", Ly = " << Ly << std::endl;
    #endif

    // Initialize u velocity (staggered at i*dx, (j+0.5)*dy)
    for (int i = 0; i < u.nx; i++) {
        for (int j = 0; j < u.ny; j++) {
            varType x = static_cast<varType>(i) * dx;
            varType y = (static_cast<varType>(j) + 0.5f) * dy;
 
            varType argSin = 2.0f * M_PI * x / Lx;
            varType argCos = 2.0f * M_PI * y / Ly;

            varType u_val = amplitude 
                          * std::sin(argSin) 
                          * std::cos(argCos);
            u.Set(i, j, u_val);
        }
    }
    
    // Initialize v velocity (staggered at (i+0.5)*dx, j*dy)
    for (int i = 0; i < u.nx; i++) {
        for (int j = 0; j < u.ny; j++) {
            const varType x = (static_cast<varType>(i) + 0.5f) * dx;
            const varType y = static_cast<varType>(j) * dy;
            
            varType argCos = 2.0f * M_PI * x / Lx;
            varType argSin = 2.0f * M_PI * y / Ly;
            const varType v_val = -amplitude 
                                * std::cos(argCos) 
                                * std::sin(argSin);
            v.Set(i, j, v_val);
        }
    }
    
}

void Fields2D::InitPotentialGradient() {
  
  varType amplitude = 1.0;
  int kx = 1; 
  int ky = 1; 

  int nxp = p.nx;
  int nyp = p.ny;

  const varType Lx = (varType)(nxp)*dx;
  const varType Ly = (varType)(nyp)*dy;

  Grid2D phi(nxp, nyp);

  for (int j = 0; j < nyp; ++j) {
    for (int i = 0; i < nxp; ++i) {
      const varType x = ((varType)i + (varType)0.5) * dx;
      const varType y = ((varType)j + (varType)0.5) * dy;
      
      varType argSinX = M_PI * (varType)kx * x / Lx;
      varType argSinY = M_PI * (varType)ky * y / Ly;
      const varType val = amplitude 
                        * std::sin(argSinX)   
                        * std::sin(argSinY);

      phi.Set(i, j, val);
    }
  }

  for (int j = 0; j < u.ny; ++j) {
    for (int i = 1; i + 1 < u.nx; ++i) {
      // u(i=0) & u(i=nx-1) already at 0
      const varType dudx = (phi.Get(i, j) - phi.Get(i - 1, j)) / dx;
      u.Set(i, j, dudx);
    }
  }

  for (int j = 1; j + 1 < v.ny; ++j) {
    for (int i = 0; i < v.nx; ++i) {
      // v(j=0) & v(j=ny-1) already at 0
      const varType dvdy = (phi.Get(i, j) - phi.Get(i, j - 1)) / dy;
      v.Set(i, j, dvdy);
    }
  }
}

void Fields2D::SolidCylinder(int cx, int cy, int r) {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      int dx = i - cx;
      int dy = j - cy;
      if (dx * dx + dy * dy <= r * r) {
        SetLabel(i, j, SOLID);
      }
    }
  }
}

void Fields2D::SolidBorders(){
  for (int i = 0; i < nx; i++) {
    SetLabel(i, 0, SOLID);
    SetLabel(i, p.ny - 1, SOLID);
  }
  for (int j = 0; j < ny; j++) {
    SetLabel(0, j, SOLID);
    SetLabel(p.nx - 1, j, SOLID);
  }
}
