#pragma once
#include "Fields.hpp"

typedef float varType;

class Project{
public:
  Project(Fields2D& fields)
    : fields(fields), nx(fields.p.nx + 1), ny(fields.p.ny + 1),
      dx(fields.dx), dy(fields.dy), dt(fields.dt),
      rhs((nx - 1) * (ny - 1), 0.0f), 
      Adiag((nx - 1) * (ny - 1), 0.0f), 
      Ax((nx - 1) * (ny - 1), 0.0f), 
      Ay((nx - 1) * (ny - 1), 0.0f) {}
 
  void MakeIncompressible();    

private:
  Fields2D& fields;
  size_t nx, ny;
  varType dx, dy, dt;
  std::vector<varType> rhs, Adiag, Ax, Ay;

  size_t idx(size_t i, size_t j) const {return (nx - 1) * j + i;}

  varType neighborPressureSum(size_t i, size_t j);

  void buildRHS();
  void buildMatrixA();
  void solveJacobi(int maxIters, varType tol);
  void updateVelocities();
};
