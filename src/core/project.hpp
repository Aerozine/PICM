#pragma once
#include "Fields.hpp"

typedef float varType;

class Project{
public:
  Project(const Fields2D& fields)
    : nx(u.nx), ny(u.ny),
      rhs((nx - 1) * (ny - 1), 0.0f), 
      Adiag((nx - 1) * (ny - 1), 0.0f), 
      Ax((nx - 1) * (ny - 1), 0.0f), 
      Ay((nx - 1) * (ny - 1), 0.0f) {}
 
  std::vector<varType> Solve();
  const std::vector<varType>& pressure() const { return p; }

private:
  const Fields2D& fields;
  size_t nx, ny;
  std::vector<varType> rhs, Adiag, Ax, Ay;

  size_t idx(size_t i, size_t j) const {return nx * j + i;}

  void buildRHS();
  void buildMatrixA();
  void solveJacobi(int maxIters, varType tol);
  varType neighborSum(size_t i, size_t j);
};
