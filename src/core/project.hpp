#pragma once
#include "Fields.hpp"

typedef float varType;

class Project{
public:
  Project(Fields2D& fields)
    : fields(fields), nx(fields.p.nx + 1), ny(fields.p.ny + 1),
      dx(fields.dx), dy(fields.dy), dt(fields.dt) {}
 
  void MakeIncompressible();    
  varType neighborPressureSum(size_t i, size_t j);
  varType neighborVelocitySum(size_t i, size_t j);
  void solveJacobi(int maxIters, float tol);
  void updateVelocities();

private:
  Fields2D& fields;
  size_t nx, ny;
  varType dx, dy, dt;
};
