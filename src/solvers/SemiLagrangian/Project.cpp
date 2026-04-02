#include "Project.hpp"

#include <iostream>

#include "SemiLagrangian.hpp"
// Pressure solve dispatch

void solvePressure(const Parameters & params,Fields2D & fields) {
  double tol = params.solver.tolerance;
  int maxIters=params.solver.maxIters;
  const double coef  = static_cast<double>(params.density) *
                       static_cast<double>(params.dx) * 
                       static_cast<double>(params.dx) /
                       static_cast<double>(params.dt);
  const double beta =  static_cast<double>(params.density) *
                       static_cast<double>(params.dx) /
                       static_cast<double>(params.dt);
  const double scale = 1.0 / coef;
  int nx =  params.nx;
  int ny =  params.ny;
  switch (params.solver.type) {
  case SolverConfig::Type::JACOBI:
    solveJacobi(fields, nx, ny, coef, maxIters, tol,beta);
    break;
  case SolverConfig::Type::GAUSS_SEIDEL:
    solveGaussSeidel(fields, nx, ny, coef, maxIters, tol,beta);
    break;
  case SolverConfig::Type::RB_GS:
    solveRedBlackGaussSeidel(fields, nx, ny, coef, maxIters, tol,beta);
    break;
  case SolverConfig::Type::MICCG0:
    solveMICCG0(fields, scale, maxIters, tol);
    break;
  default:
    std::cerr << "[SemiLagrangian] Unknown pressure solver type – aborting.\n";
    std::exit(EXIT_FAILURE);
  }
}

// Velocity correction

void updateVelocities(const Parameters & params,Fields2D & fields) {
  // Explicit pressure-gradient correction on all interior faces:
  //   u^{n+1}_{i,j} = u^*_{i,j} - (dt / (rho * dx)) * (p_{i,j} - p_{i-1,j})
  const varType coef = params.dt / (params.density * params.dx);

  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int j = 0; j < fields.u.ny; ++j) {
    for (int i = 1; i < fields.u.nx - 1; ++i) {
      if (fields.Label(i - 1, j) & Fields2D::SOLID ||
          fields.Label(i, j) & Fields2D::SOLID) {
        fields.u.Set(i, j, 0.0);
        continue;
          } else if (fields.Label(i,j)& Fields2D::BC_U) {
            continue;
          }
      fields.u.Set(i, j,
                    fields.u.Get(i, j) -
                        coef * (fields.p.Get(i, j) - fields.p.Get(i - 1, j)));
    }
  }

  OMP_PRAGMA( omp parallel for collapse(2) schedule(static))
  for (int j = 1; j < fields.v.ny - 1; ++j) {
    for (int i = 0; i < fields.v.nx; ++i) {
      if ((fields.Label(i, j - 1) & Fields2D::SOLID) ||
          (fields.Label(i, j) & Fields2D::SOLID)) {
        fields.v.Set(i, j, FIELD_USOLID);
        continue;
          } else if (fields.Label(i, j) & Fields2D::BC_V) {
            continue;
          }
      fields.v.Set(i, j,
                    fields.v.Get(i, j) -
                        coef * (fields.p.Get(i, j) - fields.p.Get(i, j - 1)));
    }
  }
}
void MakeIncompressible(const Parameters & params,Fields2D & fields){
  solvePressure(params,fields);
  updateVelocities(params, fields);
}