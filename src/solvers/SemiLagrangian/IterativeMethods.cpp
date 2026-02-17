#include "SemiLagrangian.hpp"
#include <iostream>
#include <cmath>

inline double SemiLagrangian::getUpdate(int i, int j, Real coef,
                                        double sumP, int countP)
{
  if (fields->Label(i, j) != Fields2D::FLUID) return NAN;

  sumP = 0.0;
  countP = 0;

  if (i + 1 < nx) { sumP += fields->p.Get(i + 1, j); countP++; }
  if (i - 1 >= 0) { sumP += fields->p.Get(i - 1, j); countP++; }
  if (j + 1 < ny) { sumP += fields->p.Get(i, j + 1); countP++; }
  if (j - 1 >= 0) { sumP += fields->p.Get(i, j - 1); countP++; }

  if (countP == 0) return NAN;

  double div = fields->div.Get(i, j);
  double newVal = (- coef * div + sumP) / countP;

  return newVal;
}

void SemiLagrangian::SolveJacobi(int maxIters, double tol) {
  Grid2D pNew(nx, ny);
  Real coef = density * dx * dx / dt;
  fields->Div();
  int iterations = 0;
  double sumP = 0.0;
  int countP = 0;

  for (int it = 0; it < maxIters; it++) {
    double maxDiff = 0.0;

    #pragma omp parallel for collapse(2) reduction(max:maxDiff)
    for (int i = 0; i < fields->p.nx; i++) {
      for (int j = 0; j < fields->p.ny; j++) {
        double newVal = getUpdate(i, j, coef, sumP, countP); 
        maxDiff = std::max(maxDiff, std::abs(newVal - fields->p.Get(i, j)));
        pNew.Set(i, j, newVal);
      }
    }
    
    #pragma omp parallel for collapse(2) 
    for (int i = 0; i < fields->p.nx; i++) {
      for (int j = 0; j < fields->p.ny; j++) {
        if (fields->Label(i, j) == Fields2D::FLUID) {
          fields->p.Set(i, j, pNew.Get(i, j));
        }
      }
    }

    iterations++;
    if (maxDiff < tol) {
      #ifndef NDEBUG
        std::cout << "  SolvePressure Max Diff " 
            << tol << "> maxDiff" << std::endl;
      #endif
      break;
    }
  }
  #ifndef NDEBUG
    std::cout << "  SolvePressure method converged in " 
              << iterations << " iterations" << std::endl;
  #endif
  return;
}

void SemiLagrangian::SolveGaussSeidel(int maxIters, double tol) {
  varType coef = density * dx * dx / dt;
  int iterations = 0;
  double sumP = 0.0;
  int countP = 0;
  fields->Div();

  for (int it = 0; it < maxIters; it++) {
    double maxDiff = 0.0;

    for (int i = 0; i < fields->p.nx; i++) {
      for (int j = 0; j < fields->p.ny; j++) {
        if (fields->Label(i, j) != Fields2D::FLUID) continue;

        double newVal = getUpdate(i, j, coef, sumP, countP);
        maxDiff = std::max(maxDiff, std::abs(newVal - fields->p.Get(i, j)));
        fields->p.Set(i, j, newVal);
      }
    }
    
    iterations++;
    if (maxDiff < tol) {
      #ifndef NDEBUG
        std::cout << "  SolvePressure Max Diff " 
            << tol << "> maxDiff" << std::endl;
      #endif
      break;
    }
  }
  #ifndef NDEBUG
    std::cout << "  SolvePressure method converged in " 
              << iterations << " iterations" << std::endl;
  #endif
  return;
}