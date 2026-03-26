#include "PIC.hpp"
#include "../../core/IterativeMethods.hpp"
#include <iostream>

// ---- Pressure solve dispatch ----------------------------------------------

void PIC::solvePressure(int maxIters, double tol) {
  const double coef  = static_cast<double>(density) *
                       static_cast<double>(dx) * static_cast<double>(dx) /
                       static_cast<double>(dt);
  const double scale = 1.0 / coef; // dt / (rho * dx^2)

  switch (params.solver.type) {
  case SolverConfig::Type::JACOBI:
    SolveJacobi(maxIters, tol);
    break;
  case SolverConfig::Type::GAUSS_SEIDEL:
    SolveGaussSeidel(maxIters, tol);
    break;
  case SolverConfig::Type::RB_GS:
    SolveRedBlackGaussSeidel(maxIters, tol);
    break;
  case SolverConfig::Type::MICCG0:
    solveMICCG0(*fields, scale, maxIters, tol);
    break;
  default:
    std::cerr << "[PIC] Unknown pressure solver type – aborting.\n";
    std::exit(EXIT_FAILURE);
  }
}

// ---- Velocity correction -------------------------------------------------

void PIC::updateVelocities() {
  // u^{n+1}_{i,j} = u*_{i,j} - (dt / (rho*dx)) * (p_{i,j} - p_{i-1,j})
  const varType coef = dt / (density * dx);

  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int j = 0; j < fields->u.ny; ++j) {
    for (int i = 1; i < fields->u.nx - 1; ++i) {
      if ((fields->Label(i - 1, j) & Fields2D::SOLID) ||
          (fields->Label(i,     j) & Fields2D::SOLID)) {
        fields->u.Set(i, j, fields->usolid);
        continue;
      }
      if (fields->Label(i, j) & Fields2D::BC_U)
        continue;
      fields->u.Set(i, j,
                    fields->u.Get(i, j) -
                        coef * (fields->p.Get(i, j) - fields->p.Get(i - 1, j)));
    }
  }

  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int j = 1; j < fields->v.ny - 1; ++j) {
    for (int i = 0; i < fields->v.nx; ++i) {
      if ((fields->Label(i, j - 1) & Fields2D::SOLID) ||
          (fields->Label(i, j)     & Fields2D::SOLID)) {
        fields->v.Set(i, j, fields->usolid);
        continue;
      }
      if (fields->Label(i, j) & Fields2D::BC_V)
        continue;
      fields->v.Set(i, j,
                    fields->v.Get(i, j) -
                        coef * (fields->p.Get(i, j) - fields->p.Get(i, j - 1)));
    }
  }
}

void PIC::MakeIncompressible() {
  solvePressure(params.solver.maxIters, params.solver.tolerance);
  updateVelocities();
}
