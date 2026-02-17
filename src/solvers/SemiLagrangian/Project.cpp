#include "SemiLagrangian.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

void SemiLagrangian::solvePressure(int maxIters, double tol) {
  switch (params.solver.type) {
  // TODO add preconditionning of Jacobi
  case SolverConfig::Type::JACOBI:
    SolveJacobi(maxIters, tol);
    break;
  case SolverConfig::Type::GAUSS_SEIDEL:
    SolveGaussSeidel(maxIters, tol);
    break;
  default:
    std::cerr << "Unknown pressure solver type – aborting.\n";
    exit(EXIT_FAILURE);
  }
}

void SemiLagrangian::updateVelocities() {
  varType coef = dt / (density * dx);

#pragma omp parallel for collapse(2)
  for (int i = 1; i < fields->u.nx - 1; i++) {
    for (int j = 0; j < fields->u.ny; j++) {
      if (fields->Label(i - 1, j) == Fields2D::SOLID ||
          fields->Label(i, j) == Fields2D::SOLID) {
        fields->u.Set(i, j, fields->usolid);
        continue;
      }
      varType uNew = fields->u.Get(i, j) -
                     coef * (fields->p.Get(i, j) - fields->p.Get(i - 1, j));
      fields->u.Set(i, j, uNew);
    }
  }

#pragma omp parallel for collapse(2)
  for (int i = 0; i < fields->v.nx; i++) {
    for (int j = 1; j < fields->v.ny - 1; j++) {
      if (fields->Label(i, j - 1) == Fields2D::SOLID ||
          fields->Label(i, j) == Fields2D::SOLID) {
        fields->v.Set(i, j, fields->usolid);
        continue;
      }
      varType vNew = fields->v.Get(i, j) -
                     coef * (fields->p.Get(i, j) - fields->p.Get(i, j - 1));
      fields->v.Set(i, j, vNew);
    }
  }
}

void SemiLagrangian::MakeIncompressible() {
  solvePressure(params.solver.maxIters, params.solver.tolerance);
  updateVelocities();
}
