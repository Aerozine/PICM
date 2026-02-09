#include "project.hpp"
#include <cmath>
#include <cassert>

typedef float varType;

void Project::buildRHS(){
  size_t invDx = 1 / fields.dx;

  for (size_t j = 0; j < ny - 1; j++) {
    for (size_t i = 0; i < nx - 1; i++) {
      if (fields.Label(i, j) != Fields2D::FLUID) {
        rhs[idx(i, j)] = 0.0; continue;
      }

      rhs[idx(i, j)] = -invDx * (
          (fields.u.Get(i + 1, j) - fields.u.Get(i, j)) +
          (fields.v.Get(i, j + 1) - fields.v.Get(i, j))); // eq. 5.6
    }
  }
  
  // usolid = 0 for now (adapt later)  
  for (size_t j = 0; j < ny - 1; j++) {
    for (size_t i = 0; i < nx - 1; i++) {
      
      // hardcoded version (only borders are solid) 
      /*
      if (i = 0) {
        rhs[idx(i, j)] += invDx * (fields.u.Get(i, j) - fields.usolid);
      }
      if (i = nx - 1) {
         rhs[idx(i, j)] -= invDx * (fields.u.Get(i + 1, j) - fields.usolid);
      }
      if (j = 0) {
         rhs[idx(i, j)] -= invDx * (fields.v.Get(i, j) - fields.usolid);
      }
      if (j = ny - 1) {
         rhs[idx(i, j)] += invDx * (fields.v.Get(i, j + 1) - fields.usolid);
      }
      */
      
      if (fields.Label(i, j) != Fields2D::FLUID) continue;
      
      if (i > 0 && fields.Label(i - 1, j) == Fields2D::SOLID)
         rhs[idx(i, j)] += invDx * (fields.u.Get(i, j) - fields.usolid);

      if (i + 1 < nx && fields.Label(i + 1, j) == Fields2D::SOLID)
         rhs[idx(i, j)] -= invDx * (fields.u.Get(i + 1, j) - fields.usolid);

      if (j > 0 && fields.Label(i, j - 1) == Fields2D::SOLID)
         rhs[idx(i, j)] -= invDx * (fields.v.Get(i, j) - fields.usolid);

      if (j + 1 < ny && fields.Label(i, j + 1) == Fields2D::SOLID) 
         rhs[idx(i, j)] += invDx * (fields.v.Get(i, j + 1) - fields.usolid);
    }
  }
  return;
}

void Project::buildMatrixA() {
  const double scaleA = fields.dt 
                      / (fields.density * fields.dx * fields.dx);

  for (size_t j = 0; j < ny - 1; j++) {
    for (size_t i = 0; i < nx - 1; i++) {
      if (fields.Label(i,j) != Fields2D::FLUID) continue;

      double diag = 0.0;

      // + x neighbor : (i + 1, j)
      if (i + 1 < nx && fields.Label(i + 1,j) != Fields2D::SOLID) {
          diag += scaleA;
          if (fields.Label(i + 1, j) == Fields2D::FLUID) {
            Ax[idx(i,j)] = -scaleA; // no need to do it for i - 1 (symetry)
          } 
      }
      // -x neighbor : (i - 1, j)
      if (i > 0 && fields.Label(i - 1, j) != Fields2D::SOLID) {
          diag += scaleA;
      }

      // ---- +y neighbor : (i, j+1)
      if (j + 1 < ny && fields.Label(i, j + 1) != Fields2D::SOLID) {
          diag += scaleA;
          if (fields.Label(i, j + 1) == Fields2D::FLUID) {
            Ay[idx(i, j)] = -scaleA;
          }
      }
      // ---- -y neighbor : (i, j-1)
      if (j > 0 && fields.Label(i, j - 1) != Fields2D::SOLID) {
          diag += scaleA;
      }

      Adiag[idx(i, j)] = diag;
    }
  }
  return;
}

varType Project::neighborPressureSum(size_t i, size_t j) {
  assert(i < nx - 1 && j < ny - 1); // from <cassert>
  varType sum = 0.0;

  if (i + 1 < nx - 1) sum += Ax[idx(i, j)] * fields.p.Get(i + 1, j);
  if (i > 0) sum += Ax[idx(i - 1, j)] * fields.p.Get(i - 1, j);
  if (j + 1 < ny - 1) sum += Ay[idx(i, j)] * fields.p.Get(i, j + 1);
  if (j > 0) sum += Ay[idx(i, j - 1)] * fields.p.Get(i, j - 1);
  return sum;
  // Q: neighbor pij outside domaine -> 0 contriubtion to sum ? 
}

void Project::solveJacobi(int maxIters, float tol) {
  Grid2D pNew(nx - 1, ny - 1);

  for (int it = 0; it < maxIters; it++) {

      double maxDiff = 0.0;

      for (size_t j = 0; j < ny - 1; j++) {
        for (size_t i = 0; i < nx - 1; i++) {
          if (fields.Label(i,j) != Fields2D::FLUID) continue;

          double diag = Adiag[idx(i,j)];
          if (diag == 0.0) continue;

          double sumN = neighborPressureSum(i, j);
          double newVal = (rhs[idx(i,j)] - sumN) / diag;

          maxDiff = std::max(maxDiff, std::abs(newVal - fields.p.Get(i,j)));
          pNew.Set(i,j, newVal);
        }
      }

      // swap p <- pNew (!! pNew all computed based on old values)
      for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
          if (fields.Label(i,j) == Fields2D::FLUID){  
            fields.p.Set(i,j, pNew.Get(i,j));
          }
        }
      }
      
      if (maxDiff < tol) break;
  }
  return;
}

void Project::updateVelocities() {
  
  varType coef = fields.dt / (fields.density * fields.dx);  
 
  for (size_t j = 0; j < ny - 1; j++) {
    for (size_t i = 0; i < nx; i++) {
      
      if(j + 1 == ny || j == 0) {
        fields.u.Set(i, j, fields.usolid);
        continue;
      }
      
      varType uOld = fields.u.Get(i, j); 
      varType uNew = uOld - coef * (fields.p.Get(i + 1, j) - fields.p.Get(i, j));
      fields.u.Set(i, j, uNew); 
    } 
  }

  for (size_t j = 0; j < ny; j++) {
    for (size_t i = 0; i < nx - 1; i++) {
      
      if(i + 1 == nx || i == 0) { // boundaries
        fields.v.Set(i, j, fields.usolid);
        continue;
      }
      
      // if not boundaries -> p(i, j + 1) exists ! 
      varType vOld = fields.v.Get(i, j); 
      varType vNew = vOld - coef * (fields.p.Get(i, j + 1) - fields.p.Get(i, j));
      fields.v.Set(i, j, vNew); 
    } 
  }
  return;
}

void Project::MakeIncompressible(){
  int maxIters = 10000;
  varType tol = 0.0001; 

  this -> buildRHS();
  this -> buildMatrixA();
  this -> solveJacobi(maxIters, tol); 
  this -> updateVelocities();  

  return;
}
