#include "PIC.hpp"
#include <cmath>

void PIC::Advect() {
    const varType xMax = dx * nx;
    const varType yMax = dy * ny;
    const int np = particles->size();

    keep.assign(np,1);

    fields->countAliveParticles->reset();

    OMP_PRAGMA(omp parallel for)
    for (int idx = 0; idx < np; ++idx) {
        varType x0 = particles->GetX(idx);
        varType y0 = particles->GetY(idx);
        varType u0 = particles->GetU(idx);
        varType v0 = particles->GetV(idx);

        varType xmid = x0 + static_cast<varType>(0.5) * dt * u0;
        varType ymid = y0 + static_cast<varType>(0.5) * dt * v0;

        varType umid = fields->u.interpolate<0>(xmid, ymid, dx, dy);
        varType vmid = fields->v.interpolate<1>(xmid, ymid, dx, dy);

        varType x1 = x0 + dt * umid;
        varType y1 = y0 + dt * vmid;

        if (x1 < static_cast<varType>(0) || x1 >= xMax ||
            y1 < static_cast<varType>(0) || y1 >= yMax) {
            keep[idx] = 0;
            continue;
            }

        int i1 = std::clamp(static_cast<int>(std::floor(x1 / dx)), 0, nx - 1);
        int j1 = std::clamp(static_cast<int>(std::floor(y1 / dy)), 0, ny - 1);

        if (IS_SOLID(fields->Label(i1 + 1, j1 + 1))) {
            continue; // particle stays at old position, keep[idx] remains 1
        }

        particles->SetX(idx, x1);
        particles->SetY(idx, y1);
    }

    fields->countAliveParticles->reset();

    int write = 0;
    for (int read = 0; read < np; ++read) {
        if (!keep[read]) continue;

        // compact in-place (no-op when read == write)
        if (read != write)
            particles->A[write] = particles->A[read];

        // count alive — direct struct access, no Get/Set overhead
        const varType px = particles->A[write].pos.x;
        const varType py = particles->A[write].pos.y;
        const int ci = std::clamp(static_cast<int>(std::floor(px / dx)), 0, nx - 1);
        const int cj = std::clamp(static_cast<int>(std::floor(py / dy)), 0, ny - 1);
        //@todo to be fixed
        if (!IS_SOLID(fields->Label(ci + 1, cj + 1)))
            fields->countAliveParticles->A[fields->countAliveParticles->nx * cj + ci] += 1.0f;

        ++write;
    }
    particles->A.resize(write);
}
