#include "PIC.hpp"
#include <algorithm>
#include <iostream>

static constexpr int CAPACITY_FACTOR = 3;

PIC::PIC(Parameters& params)
    : Solver(params)
{
    particles = new Particles(nx, ny, dx, dy, params.ppcx, params.ppcy);

    // Base opened the standard writers; add the PIC-only one here,
    // where the vtable is fully active and particlesWriter exists.
    if (params.write_particles)
        particlesWriter =
            std::make_unique<OutputWriter>(params.folder, "particles");

#ifndef NDEBUG
    std::cout << "Grid dimensions:\n"
              << "  p  (nx,   ny  ): " << fields->p.nx << " x " << fields->p.ny << '\n'
              << "  u  (nx+1, ny  ): " << fields->u.nx << " x " << fields->u.ny << '\n'
              << "  v  (nx,   ny+1): " << fields->v.nx << " x " << fields->v.ny << '\n';
#endif

    params.applyToFields(*fields);
    particles->InitParticleGrid(*fields);

#ifndef NDEBUG
    std::cout << "PIC initialised: " << nx << " x " << ny
              << " grid, " << params.nt << " time steps.\n\n";
#endif
}

PIC::~PIC() {
    delete particles;
    particles = nullptr;
    // fields deleted by Solver::~Solver()
}

void PIC::WriteOutput(int step) const {
    if (step % params.sampling_rate != 0)
        return;

    Solver::WriteOutput(step); // writes all standard fields

    if (params.write_particles && particlesWriter) {
        bool ok = particlesWriter->writeParticles(*particles, "particles");
        if (!ok)
            std::cerr << "[PIC] Warning: failed to write particles at step "
                      << step << '\n';
    }
}

void PIC::Step() {
    ProjectParticlesOnGrid();
    // basically a if !0 and not a inf or NAN
    if (std::isnormal(params.gravity))
        ApplyGravity();
    MakeIncompressible(params, *fields);
    ProjectGridOnParticles();
    fields->Div();
    fields->VelocityNormCenterGrid();
    Advect();
    CountAliveParticles();
    UpdateCellState();
    if (params.refill)
        RefillParticles();
}

void PIC::Run() {
    fields->Div();
    fields->VelocityNormCenterGrid();
    ProjectBCOnParticles();
    WriteOutput(0);
    RunLoop(std::max(1, params.nt / 10));
}
