#include "SemiLagrangian.hpp"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>

void SemiLagrangian::Advect() const
{
    // Create new velocity grids
    Grid2D uNew(fields->u.nx, fields->u.ny);
    Grid2D vNew(fields->v.nx, fields->v.ny);
    
    // Advect U component
    for (int i = 0; i < fields->u.nx; i++) {
        for (int j = 0; j < fields->u.ny; j++) {
            varType x, y;
            traceParticleU(i, j, x, y);
            const varType uVal = interpolateU(x, y);
            uNew.Set(i, j, uVal);
        }
    }
    
    // Advect V component
    for (int i = 0; i < fields->v.nx; i++) {
        for (int j = 0; j < fields->v.ny; j++) {
            varType x, y;
            traceParticleV(i, j, x, y);
            const varType vVal = interpolateV(x, y);
            vNew.Set(i, j, vVal);
        }
    }
    
    // Replace old velocities with new ones
    fields->u = uNew;
    fields->v = vNew;
}

void SemiLagrangian::traceParticleU(const int i, const int j,
                                    varType& x, varType& y) const
{
    // U is stored at (i*dx, (j+0.5)*dy)
    const varType x0 = static_cast<varType>(i) * dx;
    const varType y0 = (static_cast<varType>(j) + 0.5f) * dy;
    
    // RK2 backward trace - First estimate: Euler step
    varType u0, v0;
    getVelocity(x0, y0, u0, v0);
    const varType xMid = x0 - 0.5f * dt * u0;
    const varType yMid = y0 - 0.5f * dt * v0;
    
    // Second estimate: use midpoint velocity
    varType uMid, vMid;
    getVelocity(xMid, yMid, uMid, vMid);
    x = x0 - dt * uMid;
    y = y0 - dt * vMid;
    
    // Clamp to domain boundaries
    // à verifier avec domaines x,y physique et u de taille nx + 1, ny
    x = std::max(0.0f, std::min(x, static_cast<varType>(nx - 1) * dx));
    y = std::max(0.0f, std::min(y, static_cast<varType>(ny - 1) * dy));
}

void SemiLagrangian::traceParticleV(const int i, const int j,    
                                    varType& x, varType& y) const
{
    // V is stored at ((i+0.5)*dx, j*dy)
    const varType x0 = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
    const varType y0 = static_cast<varType>(j) * dy;
    
    // RK2 backward trace
    varType u0, v0;
    getVelocity(x0, y0, u0, v0);
    const varType xMid = x0 - REAL_LITERAL(0.5) * dt * u0;
    const varType yMid = y0 - REAL_LITERAL(0.5) * dt * v0;
    
    varType uMid, vMid;
    getVelocity(xMid, yMid, uMid, vMid);
    x = x0 - dt * uMid;
    y = y0 - dt * vMid;
    
    // Clamp to domain boundaries
    varType xmin = std::min(x, static_cast<varType>(nx - 1) * dx);
    varType ymin = std::min(y, static_cast<varType>(ny - 1) * dy);

    x = std::max(REAL_LITERAL(0.0), xmin);
    y = std::max(REAL_LITERAL(0.0), ymin);
}

varType SemiLagrangian::interpolateU(const varType x, const varType y) const
{
    // U is stored at positions (i*dx, (j+0.5)*dy)
    const varType i_real = x / dx;
    const varType j_real = y / dy - 0.5f;
    
    // Get integer part (cell indices)
    int i = static_cast<int>(std::floor(i_real));
    int j = static_cast<int>(std::floor(j_real));
    
    // Get fractional part for interpolation
    const varType fx = i_real - static_cast<varType>(i);
    const varType fy = j_real - static_cast<varType>(j);
    
    // Clamp indices to valid range
    i = std::max(0, std::min(i, static_cast<int>(fields->u.nx) - 2));
    j = std::max(0, std::min(j, static_cast<int>(fields->u.ny) - 2));
    
    // Bilinear interpolation
    const varType u00 = fields->u.Get(i, j);
    const varType u10 = fields->u.Get(i + 1, j);
    const varType u01 = fields->u.Get(i, j + 1);
    const varType u11 = fields->u.Get(i + 1, j + 1);

    const varType u0 = (1.0f - fx) * u00 + fx * u10;
    const varType u1 = (1.0f - fx) * u01 + fx * u11;
    
    return (1.0f - fy) * u0 + fy * u1;
}

varType SemiLagrangian::interpolateV(varType x, varType y) const
{
    // V is stored at positions ((i+0.5)*dx, j*dy)
    const varType i_real = x / dx - REAL_LITERAL(0.5);
    const varType j_real = y / dy;
    
    int i = static_cast<int>(std::floor(i_real));
    int j = static_cast<int>(std::floor(j_real));

    const varType fx = i_real - static_cast<varType>(i);
    const varType fy = j_real - static_cast<varType>(j);
    
    // Clamp indices
    i = std::max(0, std::min(i, static_cast<int>(fields->v.nx) - 2));
    j = std::max(0, std::min(j, static_cast<int>(fields->v.ny) - 2));
    
    // Bilinear interpolation
    const varType v00 = fields->v.Get(i, j);
    const varType v10 = fields->v.Get(i + 1, j);
    const varType v01 = fields->v.Get(i, j + 1);
    const varType v11 = fields->v.Get(i + 1, j + 1);

    const varType v0 = (1.0f - fx) * v00 + fx * v10;
    const varType v1 = (1.0f - fx) * v01 + fx * v11;
    
    return (1.0f - fy) * v0 + fy * v1;
}

void SemiLagrangian::getVelocity(const varType x, const varType y,
                                 varType& u, varType& v) const
{
    u = interpolateU(x, y);
    v = interpolateV(x, y);
}
