#include "Fields.hpp"
#include <cmath>

void Fields2D::Div() {

  for (size_t j = 0; j < ny - 1; j++) {
    for (size_t i = 0; i < nx - 1; i++) {
      varType dudx = (u.Get(i, j) - u.Get(i, j + 1)) / nx;
      varType dvdy = (v.Get(i, j) - v.Get(i + 1, j)) / ny;

      div.Set(i, j, dudx + dvdy);
    }
  }
  return;
}

void Fields2D::InitPotentialGradient(varType amplitude, int kx, int ky) {
  const size_t nxp = p.nx;
  const size_t nyp = p.ny;

  const varType Lx = (varType)(nxp)*dx;
  const varType Ly = (varType)(nyp)*dy;

  Grid2D phi(nxp, nyp);

  for (size_t j = 0; j < nyp; ++j) {
    for (size_t i = 0; i < nxp; ++i) {
      const varType x = ((varType)i + (varType)0.5) * dx;
      const varType y = ((varType)j + (varType)0.5) * dy;

      const varType val = amplitude * std::sin(M_PI * (varType)kx * x / Lx) *
                          std::sin(M_PI * (varType)ky * y / Ly);

      phi.Set(i, j, val);
    }
  }

  for (size_t j = 0; j < u.ny; ++j) {
    for (size_t i = 1; i + 1 < u.nx; ++i) {
      // u(i=0) & u(i=nx-1) already at 0
      const varType dudx = (phi.Get(i, j) - phi.Get(i - 1, j)) / dx;
      u.Set(i, j, dudx);
    }
  }

  for (size_t j = 1; j + 1 < v.ny; ++j) {
    for (size_t i = 0; i < v.nx; ++i) {
      // v(j=0) & v(j=ny-1) already at 0
      const varType dvdy = (phi.Get(i, j) - phi.Get(i, j - 1)) / dy;
      v.Set(i, j, dvdy);
    }
  }
}
/*
Grid2D Fields2D::VelocityNormCenterGrid() {
  Grid2D velocityNorm(nx - 1, ny - 1);
  varType uTemp, vTemp, uNorm;

  for(size_t j = 0; j < ny - 1; j++) {
    for(size_t i = 0; j < nx - 1; i++) {
      uTemp = u.Interpolate(i, j, 0);
      vTemp = v.Interpolate(i, j, 1);
      uNorm = sqrt(uTemp*uTemp + vTemp*vTemp);
      velocityNorm.Set(i, j, uNorm);
    }
  }
  return velocityNorm;
}
*/
