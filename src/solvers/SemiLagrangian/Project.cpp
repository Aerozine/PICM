#include "SemiLagrangian.hpp"
#include "../../core/IterativeMethods.hpp"
#include <iostream>

// ---- Pressure solve dispatch ----------------------------------------------

void SemiLagrangian::solvePressure(int maxIters, double tol) {
  const double coef  = static_cast<double>(density) *
                       static_cast<double>(dx) * static_cast<double>(dx) /
                       static_cast<double>(dt);
  const double scale = 1.0 / coef; // dt / (rho * dx^2)

  switch (params.solver.type) {
  case SolverConfig::Type::JACOBI:
    solveJacobi(*fields, nx, ny, coef, maxIters, tol);
    break;
  case SolverConfig::Type::GAUSS_SEIDEL:
    solveGaussSeidel(*fields, nx, ny, coef, maxIters, tol);
    break;
  case SolverConfig::Type::RB_GS:
    solveRedBlackGaussSeidel(*fields, nx, ny, coef, maxIters, tol);
    break;
  case SolverConfig::Type::MICCG0:
    solveMICCG0(*fields, scale, maxIters, tol);
    break;
  default:
    std::cerr << "[SemiLagrangian] Unknown pressure solver type – aborting.\n";
    std::exit(EXIT_FAILURE);
  }
}

// ---- Velocity correction -------------------------------------------------

void SemiLagrangian::updateVelocities() {
  // u^{n+1}_{i,j} = u*_{i,j} - (dt/(rho*dx)) * (p_{i,j} - p_{i-1,j})
  const varType coef = dt / (density * dx);

  // interior update for u
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

  // velocity at boundaries: if not fixed by BCs, copy from adjacent interior value
  for (int j = 0; j < fields->u.ny; ++j) {
    
    if (! (fields->Label(0, j) & Fields2D::BC_U))
        fields->u.Set(0, j, fields->u.Get(1, j));

    if (! (fields->Label(fields->u.nx - 1, j) & Fields2D::BC_U))
        fields->u.Set(fields->u.nx - 1, j, fields->u.Get(fields->u.nx - 2, j));
  }

  // interior update for v
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

  // velocity at boundaries: if not fixed by BCs, copy from adjacent interior value
  for (int i = 0; i < fields->v.nx; ++i) {
    
    if (! (fields->Label(i, 0) & Fields2D::BC_V))
        fields->v.Set(i, 0, fields->v.Get(i, 1));

    if (! (fields->Label(i, fields->v.ny - 1) & Fields2D::BC_V))
        fields->v.Set(i, fields->v.ny - 1, fields->v.Get(i, fields->v.ny - 2));
  }
}

void SemiLagrangian::MakeIncompressible() {
  solvePressure(params.solver.maxIters, params.solver.tolerance);
  updateVelocities();
}
