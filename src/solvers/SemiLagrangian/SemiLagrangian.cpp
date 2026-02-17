#include "SemiLagrangian.hpp"
#include <iostream>
#include <time.h>

SemiLagrangian::SemiLagrangian(const Parameters &params)
    : params(params), nx(params.nx), ny(params.ny), dx(params.dx),
      dy(params.dy), dt(params.dt), density(params.density) {
  fields = new Fields2D(nx, ny, density, dt, dx, dy);

#ifndef NDEBUG
  std::cout << "Grid dimensions:\n"
            << "  p  (nx,ny)   : " << fields->p.nx << " x " << fields->p.ny
            << '\n'
            << "  u  (nx+1,ny) : " << fields->u.nx << " x " << fields->u.ny
            << '\n'
            << "  v  (nx,ny+1) : " << fields->v.nx << " x " << fields->v.ny
            << '\n';
#endif
  // switch to select IC
  switch (params.initialCondition) {

  // WIP
  // Periodic domain: no solid walls, no cylinder.
  // Velocity is set analytically by InitTaylorGreen.
  case InitialCondition::TAYLOR_GREEN:
    fields->InitTaylorGreen(static_cast<varType>(params.taylorGreenAmplitude));
    break;

  // Custom scene
  // Apply velocity patches then solid objects from JSON.
  case InitialCondition::CUSTOM:
  default:
    for (const auto &obj : params.velocityU.objects)
      obj->applyVelocityU(*fields);
    for (const auto &obj : params.velocityV.objects)
      obj->applyVelocityV(*fields);
    for (const auto &obj : params.solid.objects)
      obj->applySolid(*fields);
    break;
  }

  InitializeOutputWriters();

#ifndef NDEBUG
  std::cout << "SemiLagrangian solver initialized: " << nx << "x" << ny
            << " grid, " << params.nt << " timesteps\n";
#endif
}

SemiLagrangian::~SemiLagrangian() { delete fields; }

void SemiLagrangian::InitializeOutputWriters() {
  if (params.write_u)
    uWriter = std::make_unique<OutputWriter>(params.folder, "u");
  if (params.write_v)
    vWriter = std::make_unique<OutputWriter>(params.folder, "v");
  if (params.write_p)
    pWriter = std::make_unique<OutputWriter>(params.folder, "p");
  if (params.write_div)
    divWriter = std::make_unique<OutputWriter>(params.folder, "div");
  if (params.write_norm_velocity)
    normVelocityWriter =
        std::make_unique<OutputWriter>(params.folder, "normVelocity");
}

void SemiLagrangian::WriteOutput(int step) const {
  if (step % params.sampling_rate != 0)
    return;

  bool ok = true;
  if (params.write_u && uWriter)
    ok &= uWriter->writeGrid2D(fields->u, "u");
  if (params.write_v && vWriter)
    ok &= vWriter->writeGrid2D(fields->v, "v");
  if (params.write_p && pWriter)
    ok &= pWriter->writeGrid2D(fields->p, "p");
  if (params.write_div && divWriter)
    ok &= divWriter->writeGrid2D(fields->div, "div");
  if (params.write_norm_velocity && normVelocityWriter)
    ok &= normVelocityWriter->writeGrid2D(fields->normVelocity, "normVelocity");

  if (!ok)
    std::cerr << "Warning: failed to write output at step " << step << '\n';
}

void SemiLagrangian::Step() {
  MakeIncompressible();
  Advect();
  fields->Div();
  fields->VelocityNormCenterGrid();
}

void SemiLagrangian::Run() {
  fields->Div();
  WriteOutput(0);

  double start = GET_TIME();

  for (int t = 1; t <= params.nt; ++t) {
    if (t % std::max(1, params.nt / 10) == 0) {
      varType max_div = 0.0f;
      for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
          max_div = std::max(max_div, std::abs(fields->div.Get(i, j)));

      std::cout << "\rStep " << t << " / " << params.nt << " ("
                << (100 * t / params.nt) << "%) "
                << "max |div| = " << max_div << std::flush;
    }

    Step();
    WriteOutput(t);
  }

  printf("\nDone: %g seconds\n", GET_TIME() - start);
}
