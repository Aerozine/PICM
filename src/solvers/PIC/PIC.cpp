#include "PIC.hpp"
#include <algorithm>
#include <iostream>

PIC::PIC(Parameters &params)
    : Solver(params), cloud(std::make_unique<Cloud2D>(params)) {
  if (params.write_particles)
    particlesWriter =
        std::make_unique<OutputWriter>(params.folder, "particles");

#ifndef NDEBUG
  std::cout << "Grid dimensions:\n"
            << "  p  (nx,   ny  ): " << fields->p.nx << " x " << fields->p.ny
            << '\n'
            << "  u  (nx+1, ny  ): " << fields->u.nx << " x " << fields->u.ny
            << '\n'
            << "  v  (nx,   ny+1): " << fields->v.nx << " x " << fields->v.ny
            << '\n';
#endif

  params.applyToFields(*fields);
  cloud->InitParticleGrid(*fields, params.ppcx, params.ppcy);

#ifndef NDEBUG
  std::cout << "PIC initialised: " << nx << " x " << ny << " grid, "
            << params.nt << " time steps, " << cloud->totalSize()
            << " particles.\n\n";
#endif
}

int PIC::computeAdvectionSubsteps() const {
  if (params.max_cfl <= REAL_EPSILON || dt <= REAL_EPSILON)
    return 1;

  varType maxCourant = varType(0);
  const varType gravitySpeed = std::abs(params.gravity) * dt;

    OMP_PRAGMA(omp parallel for collapse(2) reduction(max:maxCourant) schedule(static))
    for (int cj = 0; cj < ny; ++cj) {
      for (int ci = 0; ci < nx; ++ci) {
        const Particles &cell = (*cloud)(ci, cj);
        for (int p = 0; p < cell.size(); ++p) {
          const varType courantX = std::abs(cell.GetU(p)) / dx;
          const varType courantY = (std::abs(cell.GetV(p)) + gravitySpeed) / dy;
          maxCourant = std::max(maxCourant, std::max(courantX, courantY));
        }
      }
    }

    if (maxCourant <= REAL_EPSILON)
      return 1;

    return std::max(
        1, static_cast<int>(std::ceil((dt * maxCourant) / params.max_cfl)));
}

void PIC::WriteOutput(int step) const {
  if (step % params.sampling_rate != 0)
    return;

  if (params.write_norm_velocity)
    fields->VelocityNormCenterGrid();

  Solver::WriteOutput(step);

  if (params.write_particles && particlesWriter) {
    const bool ok = particlesWriter->writeCloud(*cloud, "particles");
    if (!ok)
      std::cerr << "[PIC] Warning: failed to write particles at step " << step
                << '\n';
  }
}

void PIC::Step() {
  const varType frameDt = dt;
  const int substeps = computeAdvectionSubsteps();
  const varType subDt = frameDt / static_cast<varType>(substeps);

  if (substeps > 1) {
    DBG_PRINTF("%s: CFL substepping %d x dt=%g",
               solverMethodName(params.solver.method).c_str(), substeps,
               static_cast<double>(subDt));
  }

  setTimeStep(subDt);
  for (int substep = 0; substep < substeps; ++substep) {
    ProjectParticlesOnGrid();
    if (params.surfaceTension && !params.particleInteraction)
      LaplacePressure();
    MakeIncompressible(params, *fields);
    ProjectGridOnParticles();
    if (params.surfaceTension && params.particleInteraction) {
      particleInteraction();
      ProjectParticlesOnGrid();
    }
    Advect();
    UpdateCellState();
    if (params.refill)
      RefillParticles();
  }
  setTimeStep(frameDt);
}

void PIC::Run() {
  fields->UpdateDivNorm();
  ProjectGridOnParticles();
  WriteOutput(0);
  RunLoop(std::max(1, params.nt / 100));
}
