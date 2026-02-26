#include "SemiLagrangian.hpp"
#include <iostream>

// Pressure solve dispatch

void SemiLagrangian::solvePressure(int maxIters, double tol) {
  switch (params.solver.type) {
  case SolverConfig::Type::JACOBI:
    SolveJacobi(maxIters, tol);
    break;
  case SolverConfig::Type::GAUSS_SEIDEL:
    SolveGaussSeidel(maxIters, tol);
    break;
  case SolverConfig::Type::RED_BLACK_GAUSS_SEIDEL:
    SolveRedBlackGaussSeidel(maxIters, tol);
    break;
  default:
    std::cerr << "[SemiLagrangian] Unknown pressure solver type – aborting.\n";
    std::exit(EXIT_FAILURE);
  }
}

// Velocity correction

void SemiLagrangian::updateVelocities() {
  // Explicit pressure-gradient correction on all interior faces:
  //   u^{n+1}_{i,j} = u^*_{i,j} - (dt / (rho * dx)) * (p_{i,j} - p_{i-1,j})
  //
  // Faces adjacent to a SOLID cell are set to usolid (no-slip wall).
  // The outermost layer of faces (i=0 and i=nx for u; j=0 and j=ny for v)
  // is left unchanged — it represents the domain boundary.

  const varType coef = dt / (density * dx);

  // u-faces: i is the fast (inner) index — contiguous in row-major storage.
#pragma omp parallel for collapse(2) schedule(static)
  for (int j = 0; j < fields->u.ny; ++j) {
    for (int i = 1; i < fields->u.nx; ++i) {
      if (fields->Label(i - 1, j) == Fields2D::SOLID ||
          fields->Label(i,     j) == Fields2D::SOLID) {
        fields->u.Set(i, j, fields->usolid);
        continue;
      }
      fields->u.Set(i, j,
                    fields->u.Get(i, j) -
                        coef * (fields->p.Get(i, j) - fields->p.Get(i - 1, j)));
    }
  }

  // v-faces: i is the fast (inner) index.
#pragma omp parallel for collapse(2) schedule(static)
  for (int j = 1; j < fields->v.ny; ++j) {
    for (int i = 0; i < fields->v.nx; ++i) {
      if (fields->Label(i, j - 1) == Fields2D::SOLID ||
          fields->Label(i, j    ) == Fields2D::SOLID) {
        fields->v.Set(i, j, fields->usolid);
        continue;
      }
      fields->v.Set(i, j,
                    fields->v.Get(i, j) -
                        coef * (fields->p.Get(i, j) - fields->p.Get(i, j - 1)));
    }
  }
}

void SemiLagrangian::MakeIncompressible() {
  solvePressure(params.solver.maxIters, params.solver.tolerance);
  updateVelocities();
}
