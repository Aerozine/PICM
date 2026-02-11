#include "project.hpp"
#include <cassert>
#include <cmath>

varType Project::neighborPressureSum(size_t i, size_t j) {
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

varType Project::neighborVelocitySum(size_t i, size_t j) {
  varType sumV = 0.0;

  if (j + 1 < ny - 1)
    sumV += fields.u.Get(i, j + 1);
  sumV -= fields.u.Get(i, j);
  if (i + 1 < nx - 1)
    sumV += fields.v.Get(i + 1, j);
  sumV -= fields.v.Get(i, j);

  return sumV;
}

void Project::solveJacobi(int maxIters, double tol) {
  Grid2D pNew(nx - 1, ny - 1);
  varType coef = fields.density * dx / dt;

  for (int it = 0; it < maxIters; it++) {

    double maxDiff = 0.0;

    for (size_t j = 0; j < ny - 1; j++) {
      for (size_t i = 0; i < nx - 1; i++) {
        if (fields.Label(i, j) != Fields2D::FLUID)
          continue;

        double sumP = neighborPressureSum(i, j);
        double sumV = neighborVelocitySum(i, j);
        double newVal = 0.25 * (coef * sumV + sumP);
        // if (i == 20 && j == 20) printf("new pressure: %f", newVal);

        maxDiff = std::max(maxDiff, std::abs(newVal - fields.p.Get(i, j)));
        pNew.Set(i, j, newVal);
      }
    }

    for (size_t j = 0; j < ny - 1; j++) {
      for (size_t i = 0; i < nx - 1; i++) {
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

  varType coef = fields.dt / (fields.density * fields.dx);

  for (size_t j = 1; j < ny - 2; j++) {
    for (size_t i = 1; i < nx - 1; i++) {

      if (j + 1 == ny || j == 0) {
        fields.u.Set(i, j, fields.usolid);
        continue;
      }

      varType uOld = fields.u.Get(i, j);
      varType uNew =
          uOld - coef * (fields.p.Get(i, j + 1) - fields.p.Get(i, j));
      fields.u.Set(i, j, uNew);
    }
  }

  for (size_t j = 1; j < ny - 1; j++) {
    for (size_t i = 1; i < nx - 2; i++) {

      if (i + 1 == nx || i == 0) { // boundaries
        fields.v.Set(i, j, fields.usolid);
        continue;
      }

      // if not boundaries -> p(i, j + 1) exists !
      varType vOld = fields.v.Get(i, j);
      varType vNew =
          vOld - coef * (fields.p.Get(i + 1, j) - fields.p.Get(i, j));
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
