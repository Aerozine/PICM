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
    if (params.write_countAliveParticles)
        countAliveParticles_writer =
            std::make_unique<OutputWriter>(params.folder, "countAliveParticles");

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
        // flatten cloud into a temporary Particles for the writer
        // @todo consider a dedicated Cloud2D writer to avoid the copy
        Particles flat(params);
        for (int ci = 0; ci < nx; ++ci)
            for (int cj = 0; cj < ny; ++cj) {
                const Particles &cell = (*cloud)(ci, cj);
                for (int p = 0; p < cell.size(); ++p)
                    flat.Add(cell.GetX(p), cell.GetY(p),
                             cell.GetU(p), cell.GetV(p), 0);
            }
        const bool ok = particlesWriter->writeParticles(flat, "particles");
        if (!ok)
            std::cerr << "[PIC] Warning: failed to write particles at step "
                      << step << '\n';
    }

    if (params.write_countAliveParticles && countAliveParticles_writer) {
        // build a temporary Grid2D of counts for the writer
        Grid2D countGrid(nx, ny);
        for (int ci = 0; ci < nx; ++ci)
            for (int cj = 0; cj < ny; ++cj)
                countGrid.Set(ci, cj, static_cast<varType>(cloud->countIn(ci, cj)));
        countAliveParticles_writer->writeGrid2D(countGrid, "countAliveParticles");
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
