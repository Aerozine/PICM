#include "SemiLagrangian.hpp"
#include <algorithm>
#include <iostream>

SemiLagrangian::SemiLagrangian(const Parameters &params)
    : params(params), nx(params.nx), ny(params.ny),
      dx(static_cast<varType>(params.dx)), dy(static_cast<varType>(params.dy)),
      dt(static_cast<varType>(params.dt)),
      density(static_cast<varType>(params.density)),
      fields(new Fields2D(nx, ny, density, dt, dx, dy)) {

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

  InitializeOutputWriters();

#ifndef NDEBUG
  std::cout << "SemiLagrangian initialised: " << nx << " x " << ny << " grid, "
            << params.nt << " time steps.\n";
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
    std::cerr << "[SemiLagrangian] Warning: failed to write output at step "
              << step << '\n';
}

void SemiLagrangian::Step() {
  MakeIncompressible(); // 1. Pressure projection: enforce div u = 0.
  Advect();             // 2. Semi-Lagrangian transport of velocity.
  fields->Div();        // } Update diagnostics used for
  fields->VelocityNormCenterGrid(); // } output and progress reporting.
}

void SemiLagrangian::Run() {
  // Compute initial diagnostics and write the t=0 snapshot.
  fields->Div();
  fields->VelocityNormCenterGrid();
  WriteOutput(0);

  const double start = GET_TIME();
  const int reportEvery = std::max(1, params.nt / 10);

  for (int t = 1; t <= params.nt; ++t) {
    // Overwrite progress line in place (~every 10 %).
    if (t % reportEvery == 0) {
      varType maxDiv = REAL_LITERAL(0.0);
      for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
          maxDiv = std::max(maxDiv, std::abs(fields->div.Get(i, j)));

      std::cout << "\rStep " << t << " / " << params.nt << " ("
                << (100 * t / params.nt) << "%) "
                << "max |div| = " << maxDiv << std::flush;
    }

    Step();
    WriteOutput(t);
  }

  std::cout << "\nDone: " << (GET_TIME() - start) << " s\n";
}
