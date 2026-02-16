#include "SemiLagrangian.hpp"
#include <cassert>
#include <cmath>
#include <stdio.h>
#include <iostream>

void SemiLagrangian::solveJacobi(int maxIters, double tol) {
  Grid2D pNew(nx, ny);
  varType coef = density * dx * dx / dt;
  fields->Div();
  int iterations = 0;

  for (int it = 0; it < maxIters; it++) {
    double maxDiff = 0.0;

    #pragma omp parallel for collapse(2) reduction(max:maxDiff)
    // update even on borders ?
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        if (fields->Label(i, j) != Fields2D::FLUID) continue;

        double sumP = 0.0;
        int countP = 0;

        if (i + 1 < nx) { sumP += fields->p.Get(i + 1, j); countP++; }
        if (i - 1 >= 0) { sumP += fields->p.Get(i - 1, j); countP++; }
        if (j + 1 < ny) { sumP += fields->p.Get(i, j + 1); countP++; }
        if (j - 1 >= 0) { sumP += fields->p.Get(i, j - 1); countP++; }

        if (countP == 0) continue;

        double div = fields->div.Get(i, j);
        double newVal = (- coef * div + sumP) / countP;

        maxDiff = std::max(maxDiff, std::abs(newVal - fields->p.Get(i, j)));
        pNew.Set(i, j, newVal);
      }
    }
    
    #pragma omp parallel for collapse(2) 
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        if (fields->Label(i, j) == Fields2D::FLUID) {
          fields->p.Set(i, j, pNew.Get(i, j));
          if (it % 499 == 0) {
            /*if (i == int(nx / 2) && j == int(ny / 2)) {
              #ifndef NDEBUG
                std::cout << "  Jacobi Iter " << it 
                          << ", p(" << i << "," << j << ") = " 
                          << fields->p.Get(i, j) << std::endl;
              #endif
            }
            */
          }
        }
      }
    }
    iterations++;
    if (maxDiff < tol) {
      /*
      #ifndef NDEBUG
        std::cout << "  Jacobi Max Diff " 
            << tol << "> maxDiff" << std::endl;
      #endif
      */
      break;
    }
  }
    /*
    #ifndef NDEBUG
        std::cout << "  Jacobi converged in " 
            << iterations << " iterations" << std::endl;
    #endif
    */
  return;
}

void SemiLagrangian::updateVelocities() {
  // to do : update velocities at borders 
  varType coef = dt / (density * dx);

  #pragma omp parallel for collapse(2) 
  for (int i = 1; i < fields->u.nx - 1; i++) {
    for (int j = 1; j < fields->u.ny - 1; j++) {

      varType uOld = fields->u.Get(i, j);
      varType uNew =
          uOld - coef * (fields->p.Get(i, j) - fields->p.Get(i - 1, j));
      fields->u.Set(i, j, uNew);
    }
  }

  #pragma omp parallel for collapse(2) 
  for (int i = 1; i < fields->v.nx - 1; i++) {
    for (int j = 1; j < fields->v.ny - 1; j++) {

      varType vOld = fields->v.Get(i, j);
      varType vNew =
          vOld - coef * (fields->p.Get(i, j) - fields->p.Get(i, j - 1));
      fields->v.Set(i, j, vNew);
    }
  }
  return;
}

void SemiLagrangian::MakeIncompressible() {
  int maxIters = 1000;
  varType tol = 1e-3;

  this->solveJacobi(maxIters, tol);
  this->updateVelocities();

  return;
}
