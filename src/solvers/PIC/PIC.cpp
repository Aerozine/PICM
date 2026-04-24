#include "PIC.hpp"
#include <algorithm>
#include <iostream>

PIC::PIC(Parameters &params)
    : Solver(params),
      cloud(std::make_unique<Cloud2D>(params))
{
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
    cloud->InitParticleGrid(*fields, params.ppcx, params.ppcy);

#ifndef NDEBUG
    std::cout << "PIC initialised: " << nx << " x " << ny << " grid, "
              << params.nt << " time steps, "
              << cloud->totalSize() << " particles.\n\n";
#endif
}

void PIC::WriteOutput(int step) const {
    if (step % params.sampling_rate != 0)
        return;

    Solver::WriteOutput(step);

    if (params.write_particles && particlesWriter) {
        const bool ok = particlesWriter->writeCloud(*cloud, "particles");
        if (!ok)
            std::cerr << "[PIC] Warning: failed to write particles at step "
                      << step << '\n';
    }
}

void PIC::Step() {
    // gravity is handled directly when setting particle speed in ProjectGridOnParticles
    ProjectParticlesOnGrid();
    MakeIncompressible(params, *fields);
    ProjectGridOnParticles();
    Advect();
    UpdateCellState();
    if (params.refill) RefillParticles();
}

void PIC::Run() {
    fields->Div();
    fields->VelocityNormCenterGrid();
    ProjectGridOnParticles();
    WriteOutput(0);

    RunLoop(std::max(1, params.nt / 100));
}
