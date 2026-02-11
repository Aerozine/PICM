#pragma once
#include "Fields.hpp"

class Project {
public:
  Project(Fields2D &fields)
      : fields(fields), nx(fields.nx), ny(fields.ny),
         dx(fields.dx), dy(fields.dy), dt(fields.dt) {}

  void MakeIncompressible();
  varType neighborPressureSum(int i, int j);
  varType neighborVelocitySum(int i, int j);
  void solveJacobi(int maxIters, double tol);
  void updateVelocities();

private:
  Fields2D &fields;
  int nx, ny;
  varType dx, dy, dt;
};
