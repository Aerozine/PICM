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
#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#define LSTART(r) _Pragma("omp parallel") { _Pragma("omp master") LIKWID_MARKER_START(r); }
#define LSTOP(r)  _Pragma("omp parallel") { _Pragma("omp master") LIKWID_MARKER_STOP(r);  }
#endif
void PIC::Step() {
  // gravity directly handle when setting particles speed
 LSTART("scatter");
  ProjectParticlesOnGrid();
  LSTOP("scatter");

  LSTART("pressure");
  MakeIncompressible(params, *fields);
  LSTOP("pressure");

  LSTART("g2p");
  ProjectGridOnParticles();
  LSTOP("g2p");

  LSTART("advect");
  Advect();
  LSTOP("advect");

  LSTART("cellstate");
  UpdateCellState();
  LSTOP("cellstate");

  if (params.refill) RefillParticles();
}

void PIC::Run() {
  fields->Div();
  fields->VelocityNormCenterGrid();
  //ProjectBCOnParticles();
  ProjectGridOnParticles();
  WriteOutput(0);
#ifdef LIKWID_PERFMON
  LIKWID_MARKER_INIT;
  OMP_PRAGMA(omp parallel)
  {
    LIKWID_MARKER_THREADINIT;
    // register all regions on all threads before the run loop
    LIKWID_MARKER_REGISTER("scatter");
    LIKWID_MARKER_REGISTER("pressure");
    LIKWID_MARKER_REGISTER("g2p");
    LIKWID_MARKER_REGISTER("advect");
    LIKWID_MARKER_REGISTER("cellstate");
  }
#endif
  RunLoop(std::max(1, params.nt / 100));
  #ifdef LIKWID_PERFMON
  LIKWID_MARKER_CLOSE;
#endif
}
