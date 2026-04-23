#include "SemiLagrangian.hpp"
#include <algorithm>
#include <iostream>

SemiLagrangian::SemiLagrangian(Parameters &params) : Solver(params) ,smokeMap(std::make_unique<Grid2D>(params.nx,params.ny)) {
#ifndef NDEBUG
  std::cout << "Grid dimensions:\n"
            << "  p  (nx,   ny  ): " << fields->p.nx << " x " << fields->p.ny
            << '\n'
            << "  u  (nx+1, ny  ): " << fields->u.nx << " x " << fields->u.ny
            << '\n'
            << "  v  (nx,   ny+1): " << fields->v.nx << " x " << fields->v.ny
            << '\n';
#endif
  params.applySmoke(std::move(*smokeMap),*fields);
  params.applyToFields(*fields);

   InitializeOutputWriters();
  if (params.write_smoke)
    smokeWriter = std::make_unique<OutputWriter>(params.folder, "smoke");
#ifndef NDEBUG
  std::cout << "SemiLagrangian initialised: " << nx << " x " << ny << " grid, "
            << params.nt << " time steps.\n";
#endif
}

void SemiLagrangian::AdvectSmoke() const {
    Grid2D smokeNew(smokeMap->nx, smokeMap->ny);

    OMP_PRAGMA(omp parallel for collapse(2))
    for (int j = 0; j < smokeMap->ny; ++j) {
        for (int i = 0; i < smokeMap->nx; ++i) {
            // if (fields->Label(i, j) & Fields2D::BC_S) {
            if (IS_BC_S(fields->Label(i + 1, j + 1))) {
                assert(std::isfinite(smokeMap->Get(i, j)));
                smokeNew.Set(i, j, smokeMap->Get(i, j));
                continue;
            }

            // Position physique du centre de la cellule (i, j)
            const varType x0 = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
            const varType y0 = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;

            // RK2 backward trace
            varType u0, v0;
            getVelocity(x0, y0, u0, v0);
            const varType xMid = x0 - REAL_LITERAL(0.5) * dt * u0;
            const varType yMid = y0 - REAL_LITERAL(0.5) * dt * v0;

            varType uMid, vMid;
            getVelocity(xMid, yMid, uMid, vMid);
            varType xDep = x0 - dt * uMid;
            varType yDep = y0 - dt * vMid;

            xDep = std::clamp(xDep, REAL_LITERAL(0.0),
                              static_cast<varType>(nx - 1) * dx);
            yDep = std::clamp(yDep, REAL_LITERAL(0.0),
                              static_cast<varType>(ny - 1) * dy);

            assert(std::isfinite(xDep));
            assert(std::isfinite(yDep));
            smokeNew.Set(i, j, smokeMap->interpolate<2>(xDep, yDep, dx, dy));
        }
    }
  // TODO : pass smokeNew as pointer of fields and avoir copy
  for (int j = 0; j < smokeMap->ny; ++j) {
    for (int i = 0; i < smokeMap->nx; ++i) {
      assert(std::isfinite(smokeNew.Get(i, j)));
      smokeMap->Set(i, j, smokeNew.Get(i, j));
    }
  }
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
  if (params.write_smoke && smokeWriter)
    smokeWriter->writeGrid2D(*smokeMap, "smoke");
  RunLoop(std::max(1, params.nt / 20));
}
