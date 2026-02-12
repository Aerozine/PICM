#include <cmath>
#include <cassert>
#include <algorithm>

void Advect(fields) {
    // new values
    Grid2D uNew(fields.u.nx, fields.u.ny);
    Grid2D vNew(fields.v.nx, fields.v.ny);
    
    // Advect for U
    for (size_t j = 0; j < fields.u.ny; ++j) {
        for (size_t i = 0; i < fields.u.nx; ++i) {
            varType x, y;
            traceParticleU(i, j, x, y);
            varType uVal = interpolateU(x, y);
            uNew.Set(i, j, uVal);
        }
    }
    // Advect for V 
    for (size_t j = 0; j < fields.v.ny; ++j) {
        for (size_t i = 0; i < fields.v.nx; ++i) {
            varType x, y;
            traceParticleV(i, j, x, y);
            varType vVal = interpolateV(x, y);
            vNew.Set(i, j, vVal);
        }
    }
    // replacing
    fields.u = uNew;
    fields.v = vNew;
}

void traceParticleU(int i,int j, varType& x, varType& y) {
    // grid value
    // allow to cast like in c but for vartype in cpp way
    varType x0 = static_cast<varType>(i) * dx;
    varType y0 = (static_cast<varType>(j) + 0.5f) * dy;
    
    // RK2 backward trace
    // First estimate: Euler step
    varType u0, v0;
    getVelocity(x0, y0, u0, v0);
    varType xMid = x0 - 0.5f * dt * u0;
    varType yMid = y0 - 0.5f * dt * v0;
    
    // Second estimate: use midpoint velocity
    varType uMid, vMid;
    getVelocity(xMid, yMid, uMid, vMid);
    x = x0 - dt * uMid;
    y = y0 - dt * vMid;
    
    // Clamp to domain boundaries
    x = std::max(0.0f, std::min(x, static_cast<varType>(nx - 1) * dx));
    y = std::max(0.0f, std::min(y, static_cast<varType>(ny - 1) * dy));
}

// same for V
void traceParticleV(size_t i, size_t j, varType& x, varType& y) {
    varType x0 = (static_cast<varType>(i) + 0.5f) * dx;
    varType y0 = static_cast<varType>(j) * dy;
    // RK2 backward trace
    varType u0, v0;
    getVelocity(x0, y0, u0, v0);
    varType xMid = x0 - 0.5f * dt * u0;
    varType yMid = y0 - 0.5f * dt * v0;
    
    varType uMid, vMid;
    getVelocity(xMid, yMid, uMid, vMid);
    x = x0 - dt * uMid;
    y = y0 - dt * vMid;
    // Clamp to domain boundaries
    x = std::max(0.0f, std::min(x, static_cast<varType>(nx - 1) * dx));
    y = std::max(0.0f, std::min(y, static_cast<varType>(ny - 1) * dy));
}

varType interpolateU(varType x, varType y) {
    // u is stored at positions (i*dx, (j+0.5)*dy)
    // Convert physical position to grid coordinates
    varType i_real = x / dx;
    varType j_real = y / dy - 0.5f;
    
    // Get integer part (cell indices)
    int i = static_cast<int>(std::floor(i_real));
    int j = static_cast<int>(std::floor(j_real));
    
    // Get fractional part for interpolation
    varType fx = i_real - static_cast<varType>(i);
    varType fy = j_real - static_cast<varType>(j);
    
    // Clamp indices to valid range
    i = std::max(0, std::min(i, static_cast<int>(fields.u.nx) - 2));
    j = std::max(0, std::min(j, static_cast<int>(fields.u.ny) - 2));
    
    // Bilinear interpolation
    varType u00 = fields.u.Get(i, j);
    varType u10 = fields.u.Get(i + 1, j);
    varType u01 = fields.u.Get(i, j + 1);
    varType u11 = fields.u.Get(i + 1, j + 1);
    
    varType u0 = (1.0f - fx) * u00 + fx * u10;
    varType u1 = (1.0f - fx) * u01 + fx * u11;
    
    return (1.0f - fy) * u0 + fy * u1;
}

varType interpolateV(varType x, varType y) {
    // v is stored at positions ((i+0.5)*dx, j*dy)
    varType i_real = x / dx - 0.5f;
    varType j_real = y / dy;
    
    int i = static_cast<int>(std::floor(i_real));
    int j = static_cast<int>(std::floor(j_real));
    
    varType fx = i_real - static_cast<varType>(i);
    varType fy = j_real - static_cast<varType>(j);
    
    // Clamp indices
    i = std::max(0, std::min(i, static_cast<int>(fields.v.nx) - 2));
    j = std::max(0, std::min(j, static_cast<int>(fields.v.ny) - 2));
    
    // Bilinear interpolation
    varType v00 = fields.v.Get(i, j);
    varType v10 = fields.v.Get(i + 1, j);
    varType v01 = fields.v.Get(i, j + 1);
    varType v11 = fields.v.Get(i + 1, j + 1);
    
    varType v0 = (1.0f - fx) * v00 + fx * v10;
    varType v1 = (1.0f - fx) * v01 + fx * v11;
    
    return (1.0f - fy) * v0 + fy * v1;
}

void getVelocity(varType x, varType y, varType& u, varType& v) {
    u = interpolateU(x, y);
    v = interpolateV(x, y);
}

