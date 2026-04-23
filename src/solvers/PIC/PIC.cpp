#include "PIC.hpp"
#include <algorithm>
#include <iostream>


PIC::PIC(Parameters &params) : Solver(params) ,particles(std::make_unique<Particles>(params)){
  //particles = new Particles(nx, ny, dx, dy, params.ppcx, params.ppcy);

  // Base opened the standard writers; add the PIC-only one here,
  // where the vtable is fully active and particlesWriter exists.
  if (params.write_particles)
    particlesWriter =
        std::make_unique<OutputWriter>(params.folder, "particles");
  if (params.write_countAliveParticles)
    countAliveParticles_writer =
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
  particles->InitParticleGrid(*fields);

#ifndef NDEBUG
  std::cout << "PIC initialised: " << nx << " x " << ny << " grid, "
            << params.nt << " time steps.\n\n";
#endif
}


void PIC::WriteOutput(int step) const {
  if (step % params.sampling_rate != 0)
    return;

  Solver::WriteOutput(step); // writes all standard fields

  if (params.write_particles && particlesWriter) {
    bool ok = particlesWriter->writeParticles(*particles, "particles");
  
    if (!ok)
      std::cerr << "[PIC] Warning: failed to write particles at step " << step
                << '\n';
  }
  if(params.write_countAliveParticles && fields->countAliveParticles)
    countAliveParticles_writer->writeGrid2D(*fields->countAliveParticles, "countAliveParticles");
}

void PIC::Step() {
  // gravity directly handle when setting particles speed
  ProjectParticlesOnGrid();
  MakeIncompressible(params, *fields);
  ProjectGridOnParticles();
  fields->UpdateDivNorm();
  Advect();
  UpdateCellState();
  if (params.refill)
    RefillParticles();
}

void PIC::Run() {
  fields->Div();
  fields->VelocityNormCenterGrid();
  //ProjectBCOnParticles();
  ProjectGridOnParticles();
  WriteOutput(0);
  RunLoop(std::max(1, params.nt / 100));
}
