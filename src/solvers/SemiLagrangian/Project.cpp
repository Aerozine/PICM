#include "SemiLagrangian.hpp"
#include <cassert>
#include <cmath>
#include <stdio.h>
#include <iostream>

void SemiLagrangian::solvePressure(int maxIters, double tol, const char* method) {

  const bool jacobi = strcmp(method, "Jacobi") == 0;
  const bool GS = strcmp(method, "Gauss-Seidel") == 0;

  if (jacobi) {
    SolveJacobi(maxIters, tol);
  } 
  else if (GS) {
    SolveGaussSeidel(maxIters, tol);
  } 
  else {
    std::cerr << "Unknown pressure solver method: " << method << std::endl;
    exit(EXIT_FAILURE);
  }
}

void SemiLagrangian::updateVelocities() {
  varType coef = dt / (density * dx);

  // update u (interior)
  #pragma omp parallel for collapse(2) 
  for (int i = 1; i < fields->u.nx - 1; i++) {
    for (int j = 0; j < fields->u.ny; j++) {
      
      if(fields->Label(i - 1, j) == Fields2D::SOLID || 
          fields->Label(i, j) == Fields2D::SOLID) {
        fields->u.Set(i, j, fields->usolid);
        continue;
      }

      varType uOld = fields->u.Get(i, j);
      varType uNew =
          uOld - coef * (fields->p.Get(i, j) - fields->p.Get(i - 1, j));
      fields->u.Set(i, j, uNew);
    }
  }

  // update v (interior)
  #pragma omp parallel for collapse(2) 
  for (int i = 0; i < fields->v.nx; i++) {
    for (int j = 1; j < fields->v.ny - 1; j++) {
      
      if(fields->Label(i - 1, j) == Fields2D::SOLID ||
         fields->Label(i, j) == Fields2D::SOLID) {
        fields->v.Set(i, j, fields->usolid);
        continue;
      }

      varType vOld = fields->v.Get(i, j);
      varType vNew =
          vOld - coef * (fields->p.Get(i, j) - fields->p.Get(i, j - 1));
      fields->v.Set(i, j, vNew);
    }
  }
}

void SemiLagrangian::MakeIncompressible(const char* method) {
  int maxIters = 1000;
  varType tol = 1e2;

  this->solvePressure(maxIters, tol, method);
  this->updateVelocities();

  return;
}
