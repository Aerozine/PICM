#include "PIC.hpp"
#include <algorithm>
#include <iostream>
static constexpr int CAPACITY_FACTOR = 3;

PIC::PIC(const Parameters &params)
    : params(params), nx(params.nx), ny(params.ny),
      dx(static_cast<varType>(params.dx)), dy(static_cast<varType>(params.dy)),
      dt(static_cast<varType>(params.dt)),
      density(static_cast<varType>(params.density)),
      fields(new Fields2D(params.nx, params.ny, params.density, params.dt, params.dx, params.dy, "PIC")),
      particles(new Particles(nx, ny, dx, dy, params.ppcx, params.ppcy,
                              CAPACITY_FACTOR * params.ppcx * params.ppcy *
                                  params.nx * params.ny)),
      deadSlots(new ParticleSlots()) {

#ifndef NDEBUG
  std::cout << "Grid dimensions:\n"
            << "  p  (nx,   ny  ): " << fields->p.nx << " x " << fields->p.ny
            << '\n'
            << "  u  (nx+1, ny  ): " << fields->u.nx << " x " << fields->u.ny
            << '\n'
            << "  v  (nx,   ny+1): " << fields->v.nx << " x " << fields->v.ny
            << '\n'
            << "  particles capacity: " << particles->capacity << '\n';
#endif

  params.applyToFields(*fields);
  particles->InitParticleGrid();

  // Kill particles inside solid cells and register all extra slots
  // (solid + overflow capacity) onto the free-list.
  InitFreeSlots();

  InitializeOutputWriters();

#ifndef NDEBUG
  std::cout << "PIC initialised: " << nx << " x " << ny << " grid, "
            << params.nt << " time steps.\n"
            << "  Free slots after init: " << deadSlots->size() << '\n';
#endif
}

PIC::~PIC() {
  delete fields;
  delete particles;
  delete deadSlots;
}

// Build the initial free-list.
//
// Two sources of free slots:
//  1. Grid slots whose home cell is solid ,these particles are killed.
//  2. Extra capacity slots (indices >= nx*ny*ppcx*ppcy) that were never
//     activated by InitParticleGrid , they are already dead.
//
// Having a large free-list from the start means the reseeder can inject
// particles at an open inlet without waiting for particles to die elsewhere.
void PIC::InitFreeSlots() {
  const int ppcx = particles->ppcx;
  const int ppcy = particles->ppcy;
  const int gridSlots = nx * ny * ppcx * ppcy;

  // 1. Kill grid slots that landed in solid cells.
  int idx = 0;
  for (int icell = 0; icell < nx; icell++) {
    for (int jcell = 0; jcell < ny; jcell++) {
      bool solid = (fields->Label(icell, jcell) & Fields2D::SOLID) != 0;
      for (int a = 0; a < ppcx; a++) {
        for (int b = 0; b < ppcy; b++) {
          if (solid) {
            particles->SetDead(idx, true);
            deadSlots->push(idx);
          }
          ++idx;
        }
      }
    }
  }

  // 2. All extra capacity slots are already dead; push them onto the list.
  for (int i = gridSlots; i < particles->capacity; i++) {
    deadSlots->push(i);
  }
}

void PIC::InitializeOutputWriters() {
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
  if (params.write_smoke)
    smokeWriter = std::make_unique<OutputWriter>(params.folder, "smoke");
  if (params.write_particles)
    particlesWriter =
        std::make_unique<OutputWriter>(params.folder, "particles");
}

void PIC::WriteOutput(int step) const {
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
  if (params.write_smoke && smokeWriter)
    ok &= smokeWriter->writeGrid2D(fields->smokeMap, "smoke");
  if (params.write_particles && particlesWriter)
    ok &= particlesWriter->writeParticles(*particles, "particles");
  if (!ok)
    std::cerr << "[PIC] Warning: failed to write output at step " << step
              << '\n';
}

void PIC::Step() {
  ProjectParticlesOnGrid("hat");
  MakeIncompressible(params,*fields);
  fields->Div();
  fields->VelocityNormCenterGrid();
  ProjectGridOnParticles();
  AdvectParticles();
  RefillParticles();
}

void PIC::Run() {
  fields->Div();
  fields->VelocityNormCenterGrid();
  ProjectGridOnParticles();
  WriteOutput(0);

  const double start = GET_TIME();
  const int reportEvery = std::max(1, params.nt / 10);

  for (int t = 1; t <= params.nt; ++t) {
    if (t % reportEvery == 0) {
      varType maxDiv = REAL_LITERAL(0.0);
      for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
          maxDiv = std::max(maxDiv, std::abs(fields->div.Get(i, j)));

      std::cout << "\rStep " << t << " / " << params.nt << " ("
                << (100 * t / params.nt) << "%) "
                << "max |div| = " << maxDiv
                << "  free slots: " << deadSlots->size() << std::flush;
    }

    Step();
    WriteOutput(t);
  }

  std::cout << "\nDone: " << (GET_TIME() - start) << " s\n";
}
