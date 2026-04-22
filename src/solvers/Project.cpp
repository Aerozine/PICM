#include "Solver.hpp"
#include <iostream>
// Pressure solve dispatch

void Solver::solvePressure(const Parameters &p, Fields2D &f) {
  double tol = p.solver.tolerance;
  int maxIters = p.solver.maxIters;
  const double coef = static_cast<double>(p.density) *
                      static_cast<double>(p.dx) * static_cast<double>(p.dx) /
                      static_cast<double>(p.dt);
  const double beta = static_cast<double>(p.density) *
                      static_cast<double>(p.dx) / static_cast<double>(p.dt);
  const double scale = 1.0 / coef;
  int nx = p.nx;
  int ny = p.ny;
  switch (p.solver.type) {
    // @todo fix JACOBI , GS , RM MICCG0 and do CG
  case SolverConfig::Type::JACOBI:
    solveJacobi(f, nx, ny, coef, maxIters, tol, beta);
    break;
  case SolverConfig::Type::GAUSS_SEIDEL:
    solveGaussSeidel(f, nx, ny, coef, maxIters, tol, beta);
    break;
    case SolverConfig::Type::RB_GS:
#ifdef USE_CUDA
      //std::cout << "[Solver] RB_GS → GPU (CUDA)\n";
      solveRedBlackGaussSeidel_GPU(f, nx, ny, coef, maxIters, tol, beta);
#else
      solveRedBlackGaussSeidel(f, nx, ny, coef, maxIters, tol, beta);
#endif
      break;
  case SolverConfig::Type::MICCG0:
    solveMICCG0(f, scale, maxIters, tol);
    break;
    case SolverConfig::Type::CG:
#ifdef USE_CUDA
      //solveCG_GPU(f, coef, beta, maxIters, tol);
#else
      solveCG(f, coef, beta, maxIters, tol);
#endif
      break;
    default:
    std::cerr << "[SemiLagrangian] Unknown pressure solver type – aborting.\n";
    std::exit(EXIT_FAILURE);
  }
}

// Velocity correction

void Solver::updateVelocities(const Parameters &params, Fields2D &fields) {
  // Explicit pressure-gradient correction on all interior faces:
  //   u^{n+1}_{i,j} = u^*_{i,j} - (dt / (rho * dx)) * (p_{i,j} - p_{i-1,j})
  const varType coef = params.dt / (params.density * params.dx);

  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int j = 0; j < fields.u.ny; ++j) {
    for (int i = 0; i < fields.u.nx; ++i) {
      if (IS_BC_U(fields.Label(i, j + 1))) {
        continue;
      } else if (IS_SOLID(fields.Label(i + 1, j + 1)) ||
                 IS_SOLID(fields.Label(i, j + 1))) {
        fields.u.Set(i, j, 0.0);
        continue;
      }
      fields.u.Set(i, j,
                   fields.u.Get(i, j) - coef * (fields.p.Get(i + 1, j + 1) -
                                                fields.p.Get(i, j + 1)));
    }
  }

  OMP_PRAGMA( omp parallel for collapse(2) schedule(static))
  for (int j = 0; j < fields.v.ny; ++j) {
    for (int i = 0; i < fields.v.nx; ++i) {
      // bc is applied on p1,0 for v0,0
      if (IS_BC_V(fields.Label(i + 1, j))) {
        continue;
      } else if (IS_SOLID(fields.Label(i + 1, j + 1)) ||
                 IS_SOLID(fields.Label(i + 1, j))) {
        fields.v.Set(i, j, 0.0);
        continue;
      }
      assert(std::isfinite(
          fields.v.Get(i, j) -
          coef * (fields.p.Get(i + 1, j + 1) - fields.p.Get(i + 1, j))));
      fields.v.Set(i, j,
                   fields.v.Get(i, j) - coef * (fields.p.Get(i + 1, j + 1) -
                                                fields.p.Get(i + 1, j)));
    }
  }
}

void Solver::MakeIncompressible(const Parameters &p, Fields2D &f) {
  solvePressure(p, f);
  updateVelocities(p, f);
}
