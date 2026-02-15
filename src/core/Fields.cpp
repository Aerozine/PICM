#include "Fields.hpp"
#include <ostream>
#include <iostream>
#include <cmath>

void Fields2D::Div() {
  for (int j = 0; j < ny - 1; j++) {
    for (int i = 0; i < nx - 1; i++) {
      const varType dudx = (u.Get(i + 1, j) - u.Get(i, j)) / dx;
      const varType dvdy = (v.Get(i, j + 1) - v.Get(i, j)) / dy;
      div.Set(i, j, dudx + dvdy);
    }
  }
}

void Fields2D::InitPotentialGradient(varType amplitude, int kx, int ky) {
  const int nxp = p.nx;
  const int nyp = p.ny;

  const varType Lx = static_cast<varType>(nxp)*dx;
  const varType Ly = static_cast<varType>(nyp)*dy;

  Grid2D phi(nxp, nyp);

  for (int j = 0; j < nyp; ++j) {
    for (int i = 0; i < nxp; ++i) {
      const varType x = (static_cast<varType>(i) + static_cast<varType>(0.5)) * dx;
      const varType y = (static_cast<varType>(j) + static_cast<varType>(0.5)) * dy;

      const varType val = amplitude * std::sin(M_PI * static_cast<varType>(kx) * x / Lx) *
                          std::sin(M_PI * static_cast<varType>(ky) * y / Ly);

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
// FANCY CONFIG 
// simply seting the value of u as in the center
void Fields2D::InitTaylorGreen(const varType amplitude) {
    const varType Lx = static_cast<varType>(nx - 1) * dx;
    const varType Ly = static_cast<varType>(ny - 1) * dy;
#ifdef NDEBUG
    std::cout << "Initializing Taylor-Green vortex with amplitude " << amplitude << std::endl;
    std::cout << "Domain: Lx = " << Lx << ", Ly = " << Ly << std::endl;
#endif

    // Initialize u velocity (staggered at i*dx, (j+0.5)*dy)
    for (int j = 0; j < u.ny; ++j) {
        for (int i = 0; i < u.nx; ++i) {
            varType x = static_cast<varType>(i) * dx;
            varType y = (static_cast<varType>(j) + 0.5f) * dy;

            varType u_val = amplitude * std::sin(2.0f * M_PI * x / Lx) *
              std::cos(2.0f * M_PI * y / Ly);
            u.Set(i, j, u_val);
        }
    }
    
    // Initialize v velocity (staggered at (i+0.5)*dx, j*dy)
    for (int j = 0; j < v.ny; ++j) {
        for (int i = 0; i < v.nx; ++i) {
            const varType x = (static_cast<varType>(i) + 0.5f) * dx;
            const varType y = static_cast<varType>(j) * dy;

            const varType v_val = -amplitude * std::cos(2.0f * M_PI * x / Lx) *
              std::sin(2.0f * M_PI * y / Ly);
            v.Set(i, j, v_val);
        }
    }
    
}

