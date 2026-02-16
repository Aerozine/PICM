#include "SemiLagrangian.hpp"
#include <iostream>

SemiLagrangian::SemiLagrangian(const Parameters &params)
    : params(params),
      nx(params.nx),
      ny(params.ny),
      dx(params.dx),
      dy(params.dy),
      dt(params.dt),
      density(params.density)
{
    // Create fields
    fields = new Fields2D(nx, ny, density, dt, dx, dy);

#ifndef NDEBUG
    std::cout << "Grid dimensions:" << std::endl;
    std::cout << "  Pressure (nx, ny): " << fields->p.nx 
                << " x " << fields->p.ny << std::endl;
    std::cout << "  U velocity (nx + 1, ny): " << fields->u.nx 
                << " x " << fields->u.ny << std::endl;
    std::cout << "  V velocity (nx, ny + 1): " << fields->v.nx 
                << " x " << fields->v.ny << std::endl;
#endif

    // Initialize Taylor-Green vortex
    // fields->InitTaylorGreen(REAL_LITERAL(1.0));
    fields -> u.InitRectangle(1.0);

    // Initialize output writers
    InitializeOutputWriters();

#ifndef NDEBUG
    std::cout << "SemiLagrangian solver initialized: "
              << nx << "x" << ny << " grid, "
              << params.nt << " timesteps" << std::endl;
#endif
}

SemiLagrangian::~SemiLagrangian()
{
    delete fields;
}

void SemiLagrangian::InitializeOutputWriters()
{
    if (params.write_u)
    {
        uWriter = std::make_unique<OutputWriter>(params.folder, "u");
    }
    if (params.write_v)
    {
        vWriter = std::make_unique<OutputWriter>(params.folder, "v");
    }
    if (params.write_p)
    {
        pWriter = std::make_unique<OutputWriter>(params.folder, "p");
    }
    if (params.write_div)
    {
        divWriter = std::make_unique<OutputWriter>(params.folder, "div");
    }
}

void SemiLagrangian::WriteOutput(int step) const
{
    // Only write on sampling intervals
    if (step % params.sampling_rate != 0)
    {
        return;
    }

    bool success = true;

    if (params.write_u && uWriter)
    {
        success &= uWriter->writeGrid2D(fields->u, "u");
    }
    if (params.write_v && vWriter)
    {
        success &= vWriter->writeGrid2D(fields->v, "v");
    }
    if (params.write_p && pWriter)
    {
        success &= pWriter->writeGrid2D(fields->p, "p");
    }
    if (params.write_div && divWriter)
    {
        success &= divWriter->writeGrid2D(fields->div, "div");
    }

    if (!success)
    {
        std::cerr << "Warning: Failed to write output at step " 
                    << step << std::endl;
    }
}

void SemiLagrangian::Step()
{   
    // Vincent: jamais appliquer advect is pas un divergence free field

    // 1. Make velocity field incompressible
    MakeIncompressible();
    
    // 2. Advect velocity field
    Advect();

    // 3. Update divergence field for diagnostics
    fields->Div();
}

void SemiLagrangian::Run()
{

#ifndef NDEBUG
    std::cout << "Starting Taylor-Green vortex simulation..." << std::endl;
#endif
    // Write initial condition
    fields->Div();

    std::cout << "Writed initial cond" << std::endl;
    // Check initial divergence
    varType max_div = 0.0f;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            max_div = std::max(max_div, std::abs(fields->div.Get(i, j)));
        }
    }

    WriteOutput(0);

    // Main time-stepping loop
    for (int t = 1; t <= params.nt; ++t)
    {
        // Progress indicator with divergence check
        if (t % std::max(1, params.nt / 10) == 0)
        {
            // Check divergence
            varType max_div_step = 0.0f;
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    max_div_step = std::max(max_div_step
                                , std::abs(fields->div.Get(i, j)));
                }
            }
            // fancy printing like in HPC
            std::cout << "\rStep " << t << " / " << params.nt
                      << " (" << (100 * t / params.nt) << "%) "
                      << "max |div| = " << max_div_step
                      << std::flush;
        }

        // Perform one time step
        Step();

        // Write output
        WriteOutput(t);
    }
}
