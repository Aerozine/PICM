#include "project.hpp"
#include <cassert>
#include <cmath>
#include <stdio.h>

varType Project::neighborPressureSum(int i, int j) {
  varType sumP = 0.0;

  if (i + 1 < nx - 1)
    sumP += fields.p.Get(i + 1, j);
  if (i > 0)
    sumP += fields.p.Get(i - 1, j);
  if (j + 1 < ny - 1)
    sumP += fields.p.Get(i, j + 1);
  if (j > 0)
    sumP += fields.p.Get(i, j - 1);

  return sumP;
}

varType Project::neighborVelocitySum(int i, int j) {
  varType sumV = 0.0;

  sumV -= fields.u.Get(i, j);
  sumV -= fields.v.Get(i, j);

  if (i + 1 < nx - 1)
    sumV += fields.u.Get(i + 1, j);
  if (j + 1 < ny - 1)
    sumV += fields.v.Get(i, j + 1);

  return sumV;
}

void Project::solveJacobi(int maxIters, double tol) {
  Grid2D pNew(nx, ny);
  varType coef = fields.density * dx / dt;

  for (int it = 0; it < maxIters; it++) {

    double maxDiff = 0.0;
    
    // update even on borders ?
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        if (fields.Label(i, j) != Fields2D::FLUID) continue;

        double sumP = neighborPressureSum(i, j);
        double sumV = neighborVelocitySum(i, j);
        double newVal = 0.25 * (- coef * sumV + sumP);

        if (i == nx / 2 && j == ny / 2) printf("new pressure: %f\n", newVal);

        maxDiff = std::max(maxDiff, std::abs(newVal - fields.p.Get(i, j)));
        pNew.Set(i, j, newVal);
      }
    }

    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        if (fields.Label(i, j) == Fields2D::FLUID) {
          fields.p.Set(i, j, pNew.Get(i, j));
        }
      }
    }

    if (maxDiff < tol)
      break;
  }
  return;
}

void Project::updateVelocities() {
  // to do : update velocities at borders 
  varType coef = fields.dt / (fields.density * fields.dx);

  for (int i = 1; i < nx; i++) {
    for (int j = 1; j < nx - 1; j++) {

      varType uOld = fields.u.Get(i, j);
      varType uNew =
          uOld - coef * (fields.p.Get(i, j) - fields.p.Get(i - 1, j));
      fields.u.Set(i, j, uNew);
    }
  }

  for (int i = 1; i < nx - 1; i++) {
    for (int j = 1; j < ny; j++) {

      varType vOld = fields.v.Get(i, j);
      varType vNew =
          vOld - coef * (fields.p.Get(i, j) - fields.p.Get(i, j - 1));
      fields.v.Set(i, j, vNew);
    }
  }
  return;
}

void Project::MakeIncompressible() {
  int maxIters = 10;
  varType tol = 0.0001;

  this->solveJacobi(maxIters, tol);
  this->updateVelocities();

  return;
}
