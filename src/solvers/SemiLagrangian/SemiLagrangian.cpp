#include "SemiLagrangian.hpp"
#include <algorithm>
#include <iostream>

SemiLagrangian::SemiLagrangian(Parameters &params) : Solver(params) {
/*params(params), nx(params.nx), ny(params.ny),
      dx(static_cast<varType>(params.dx)), dy(static_cast<varType>(params.dy)),
      dt(static_cast<varType>(params.dt)),
      density(static_cast<varType>(params.density)),
      fields(
          new Fields2D(nx, ny, density, dt, dx, dy, "SL", params.freeSurface)) {
*/
#ifndef NDEBUG
  std::cout << "Grid dimensions:\n"
            << "  p  (nx,   ny  ): " << fields->p.nx << " x " << fields->p.ny
            << '\n'
            << "  u  (nx+1, ny  ): " << fields->u.nx << " x " << fields->u.ny
            << '\n'
            << "  v  (nx,   ny+1): " << fields->v.nx << " x " << fields->v.ny
            << '\n';
#endif

  // Apply initial conditions from the JSON config (velocity patches, solid
  // geometry). SceneObject instances are created and destroyed inside here.
  params.applyToFields(*fields);

  // InitializeOutputWriters();

#ifndef NDEBUG
  std::cout << "SemiLagrangian initialised: " << nx << " x " << ny << " grid, "
            << params.nt << " time steps.\n";
#endif
}

void SemiLagrangian::Step() {
  MakeIncompressible(params, *fields);
  fields->Div();
  fields->VelocityNormCenterGrid();
  Advect();
  AdvectSmoke();
  DBG_PRINTF("step");
}

void SemiLagrangian::Run() {
  fields->Div();
  fields->VelocityNormCenterGrid();
  WriteOutput(0);
  RunLoop(std::max(1, params.nt / 20));
}
