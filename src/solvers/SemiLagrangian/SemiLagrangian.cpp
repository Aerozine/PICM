#include "SemiLagrangian.hpp"
#include <iostream>
#include <time.h>

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
    
    if (params.velocity_u) {
        fields -> u.InitRectangle(20.0);
    }
    if (params.velocity_v) {
        fields -> v.InitRectangle(20.0);
    }
    if (params.solid_cylinder) {
        fields -> SolidCylinder(nx/2 + 20, ny/2 + 20, 10);
    }
    if (params.solid_borders) {
        fields -> SolidBorders();
    }

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
    if (params.write_norm_velocity)
    {
        normVelocityWriter = std::make_unique<OutputWriter>(params.folder, "normVelocity");
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
    if (params.write_norm_velocity && normVelocityWriter)
    {
        success &= normVelocityWriter->writeGrid2D(fields->normVelocity, "normVelocity");
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
    //const char* method = "Jacobi";
    const char* method = "Gauss-Seidel";
    MakeIncompressible(method);
    
    // 2. Advect velocity field
    Advect();

    // 3. Update divergence field for diagnostics
    fields->Div();

    // 4. Compute norm of velocity at grid center
    fields->VelocityNormCenterGrid();
}

void SemiLagrangian::Run()
{
    // Write initial condition
    fields->Div();

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
    
    double start = GET_TIME();

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
    
    double time = GET_TIME() - start;
    printf("\nDone: %g seconds\n", time);
}
