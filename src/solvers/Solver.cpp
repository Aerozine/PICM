#include "Solver.hpp"

#include <iostream>

Solver::Solver(Parameters &params)
    : params(params), nx(params.nx), ny(params.ny), dx(params.dx),
      dy(params.dy), dt(params.dt), density(params.density),
      fields(new Fields2D(nx, ny, density, dt, dx, dy, params.freeSurface)) {
  InitializeOutputWriters();
}

void Solver::setTimeStep(varType newDt) noexcept {
  dt = newDt;
  params.dt = newDt;
  fields->dt = newDt;
}

void Solver::RunLoop(int reportEvery) {
  const varType start = GET_TIME();

  for (int t = 1; t <= params.nt; ++t) {
    if (t % reportEvery == 0) {
      varType maxDiv = REAL_LITERAL(0.0);
      int maxi = -1;
      int maxj = -1;
      OMP_PRAGMA(omp parallel for reduction(max:maxDiv))
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
          varType a = std::abs(fields->div.Get(i, j));
          if (a > maxDiv) {
            maxDiv = a;
            maxi = i;
            maxj = j;
          }
        }
      std::cout << "\rStep " << t << " / " << params.nt << " ("
                << (100 * t / params.nt) << "%) " << "max |div| = " << maxDiv
                << " reached at (" << maxi << "," << maxj << ")" << std::flush;
    }
    Step();
    WriteOutput(t);
  }

  std::cout << "\nDone: " << (GET_TIME() - start) << " s\n";
}

void Solver::InitializeOutputWriters() {
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
  labelWriter = std::make_unique<OutputWriter>(params.folder, "label");
}
void Solver::WriteOutput(int step) const {
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

  ok &= labelWriter->writeLabels(fields->Labels, fields->nx + 2, fields->ny + 2,
                                 "label");
  if (!ok)
    std::cerr << "[SemiLagrangian] Warning: failed to write output at step "
              << step << '\n';
}
