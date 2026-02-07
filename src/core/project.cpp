#include "project.hpp"

typedef float varType;

void Project::buildRHS(){
  size_t invDx = 1 / grid.dx;

  for (int j = 0; j < ny - 1; j++) {
    for (int i = 0; i < nx - 1; i++) {
      if (fields.Label(i, j) != FLUID) {rhs(i, j) = 0.0; continue;}

      rhs(i, j) = -invDx * (
          (u.Get(i + 1, j) - u.Get(i, j)) +
          (v.Get(i,j+1,k) - v.Get(i,j,k))); // eq. 5.6
    }
  }
  
  // usolid = 0 for now (adapt later)  
  for (int j = 0; j < ny - 1; j++) {
    for (int i = 0; i < nx - 1; i++) {
      if (fields.Label(i, j) != FLUID) continue;

      if (fields.Label(i - 1, j) == SOLID)
         rhs(i, j) += invDx * (u(i, j) - fields.usolid);

      if (fields.Label(i + 1, j) == SOLID)
         rhs(i, j) -= invDx * (u(i + 1, j) - fields.usolid);

      if (fields.Label(i, j - 1) == SOLID)
         rhs(i, j) -= invDx * (v(i, j) - fields.usolid);

      if (fields.Label(i, j + 1) == SOLID) 
         rhs(i, j) += invDx * (v(i, j + 1) - fields.usolid);
    }
  }
  return;
}

void Project::buildMatrixA() {
  const double scaleA = dt / (density * grid.dx * grid.dx);

  for (int j = 0; j < ny - 1; ++j) {
    for (int i = 0; i < nx - 1; ++i) {
      if (fields.Label(i,j) != FLUID) continue;

      double diag = 0.0;

      // + x neighbor : (i + 1, j)
      if (i + 1 < nx && fields.Label(i + 1,j) != SOLID) {
          diag += scaleA;
          if (fields.Label(i + 1, j) == FLUID) Axp(i,j) = -scaleA; // no need to do it
                                                                  // for i - 1 (symetry)
      }
      // -x neighbor : (i - 1, j)
      if (i - 1 >= 0 && fields.Label(i - 1, j) != SOLID) {
          diag += scaleA;
      }

      // ---- +y neighbor : (i, j+1)
      if (j + 1 < ny && fields.Label(i, j + 1) != SOLID) {
          diag += scaleA;
          if (fields.Label(i, j + 1) == FLUID) Ay(i, j) = -scaleA;
      }
      // ---- -y neighbor : (i, j-1)
      if (j - 1 >= 0 && fields.Label(i, j - 1) != SOLID) {
          diag += scaleA;
      }

      Adiag(idx(i, j)) = diag;
    }
  }
}

varType Project::neighborPressureSum(size_t i, size_t j) {
    varType sum = 0.0;
    if (i + 1 < nx) sum += Ax(idx(i, j)) * p(idx(i + 1, j));
    if (i - 1 >= 0) sum += Ax(idx(i - 1, j)) * p(idx(i - 1, j));
    if (j + 1 < ny) sum += Ay(idx(i, j)) * p(idx(i, j + 1));
    if (j - 1 >= 0) sum += Ay(idx(i, j - 1)) * p(idx(i, j - 1));
    return s;
}

void Project::solveJacobi(int maxIters, float tol) {
  Grid2D pNew(nx - 1, ny - 1);

  for (int it = 0; it < maxIters; it++) {

      double maxDiff = 0.0;

      for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
          if (label(i,j) != FLUID) continue;

          double diag = Adiag(i,j);
          if (diag == 0.0) continue;

          double sumN = neighborSum(i, j);
          double newVal = (rhs(i,j) - sumN) / diag;

          maxDiff = std::max(maxDiff, std::abs(newVal - p(i,j)));
          pNew(i,j) = newVal;
        }
      }

      // swap p <- pNew (!! pNew all computed based on old values)
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          if (label(i,j) == FLUID) p(i,j) = pNew(i,j);
        }
      }
      
      if (maxUpdate < tol) break;
  }
}

