#include "Solver.hpp"
#include <iostream>
// Pressure solve dispatch
// select the solver based on the solverconfig
inline void solvePressure(const Parameters &p, Fields2D &f) {
  varType tol = p.solver.tolerance;
  int maxIters = p.solver.maxIters;
  const varType coef = p.density * p.dx * p.dx / p.dt;
  const varType beta = p.density * p.dx / p.dt;
  int nx = p.nx;
  int ny = p.ny;
  switch (p.solver.type) {
  case SolverConfig::Type::JACOBI:
    solveJacobi(f, nx, ny, coef, maxIters, tol, beta);
    break;
  case SolverConfig::Type::GAUSS_SEIDEL:
    solveGaussSeidel(f, nx, ny, coef, maxIters, tol, beta);
    break;
  case SolverConfig::Type::RB_GS:
#ifdef USE_CUDA
    solveRedBlackGaussSeidel_GPU(f, nx, ny, coef, maxIters, tol, beta);
#else
    solveRedBlackGaussSeidel(f, nx, ny, coef, maxIters, tol, beta);
#endif
    break;
  case SolverConfig::Type::MICCG0:
    solveMICCG0(f, coef, maxIters, tol);
    break;
  case SolverConfig::Type::CG:
#ifdef USE_CUDA
    solveCG_GPU(f, coef, beta, maxIters, tol);
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

inline void updateVelocities(const Parameters &params, Fields2D &fields) {
  // Explicit pressure-gradient correction on all interior faces:
  //   u^{n+1}_{i,j} = u^*_{i,j} - (dt / (rho * dx)) * (p_{i,j} - p_{i-1,j})
  const varType coef = params.dt / (params.density * params.dx);

  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int j = 0; j < fields.u.ny; ++j) {
    for (int i = 0; i < fields.u.nx; ++i) {

      labeltype left = fields.Label(i, j + 1);
      labeltype right = fields.Label(i + 1, j + 1);

      if (IS_BC_U(left)) {
        continue;
      } else if (IS_SOLID(right) || IS_SOLID(left)) {
        fields.u.Set(i, j, 0.0);
        continue;
      }

      varType pLeft = fields.p.Get(i, j + 1);
      varType pRight = fields.p.Get(i + 1, j + 1);

      if (params.surfaceTension) {
        if (IS_AIR(left) && IS_FLUID(right)) {
          pLeft = fields.interface_u->Get(i, j);
        } else if (IS_FLUID(left) && IS_AIR(right)) {
          pRight = fields.interface_u->Get(i, j);
        }
      }
      fields.u.Set(i, j, fields.u.Get(i, j) - coef * (pRight - pLeft));
    }
  }

  OMP_PRAGMA( omp parallel for collapse(2) schedule(static))
  for (int j = 0; j < fields.v.ny; ++j) {
    for (int i = 0; i < fields.v.nx; ++i) {

      labeltype down = fields.Label(i + 1, j);
      labeltype up = fields.Label(i + 1, j + 1);

      if (IS_BC_V(down)) {
        continue;
      } else if (IS_SOLID(up) || IS_SOLID(down)) {
        fields.v.Set(i, j, 0.0);
        continue;
      }

      varType pDown = fields.p.Get(i + 1, j);
      varType pUp = fields.p.Get(i + 1, j + 1);

      if (params.surfaceTension) {
        if (IS_AIR(down) && IS_FLUID(up)) {
          pDown = fields.interface_v->Get(i, j);
        } else if (IS_FLUID(down) && IS_AIR(up)) {
          pUp = fields.interface_v->Get(i, j);
        }
      }
      fields.v.Set(i, j, fields.v.Get(i, j) - coef * (pUp - pDown));
    }
  }
}

void Solver::MakeIncompressible(const Parameters &p, Fields2D &f) {
  solvePressure(p, f);
  updateVelocities(p, f);
  f.Div();
}
