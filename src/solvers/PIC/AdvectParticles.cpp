#include "PIC.hpp"
#include <cmath>

void PIC::Advect() {
    const varType xMax = dx * nx;
    const varType yMax = dy * ny;
    const int np = particles->size();
    keep.assign(np,1);

    OMP_PRAGMA(omp parallel for schedule(static))
    for (int idx = 0; idx < np; ++idx) {
        const varType x0 = particles->GetX(idx);
        const varType y0 = particles->GetY(idx);
        const varType u0 = particles->GetU(idx);
        const varType v0 = particles->GetV(idx);

        const varType xmid = x0 + varType(0.5) * dt * u0;
        const varType ymid = y0 + varType(0.5) * dt * v0;

        const varType umid = fields->u.interpolate<0>(xmid, ymid, dx, dy);
        const varType vmid = fields->v.interpolate<1>(xmid, ymid, dx, dy);

        const varType x1 = x0 + dt * umid;
        const varType y1 = y0 + dt * vmid;

        if (x1 < varType(0) || x1 >= xMax || y1 < varType(0) || y1 >= yMax) {
            keep[idx] = 0;
            continue;
        }
        const int i1 = std::clamp(static_cast<int>(std::floor(x1 / dx)), 0, nx - 1);
        const int j1 = std::clamp(static_cast<int>(std::floor(y1 / dy)), 0, ny - 1);

        if (IS_SOLID(fields->Label(i1 + 1, j1 + 1)))
            continue; // stays at old position, keep[idx] = 1

        particles->SetX(idx, x1);
        particles->SetY(idx, y1);
    }

    // ── Compact dead particles (updates all SoA arrays atomically) ───────────
    particles->Compact(keep);

    // ── Count alive particles per cell (now trivially parallel) ─────────────
    fields->countAliveParticles->reset();
    const int np_alive = particles->size();

    // Using atomics here is fine: countAliveParticles is a small integer grid
    // and this is just an increment, not a full scatter.
    OMP_PRAGMA(omp parallel for schedule(static))
    for (int idx = 0; idx < np_alive; ++idx) {
        const int ci = std::clamp(
            static_cast<int>(std::floor(particles->GetX(idx) / dx)), 0, nx - 1);
        const int cj = std::clamp(
            static_cast<int>(std::floor(particles->GetY(idx) / dy)), 0, ny - 1);
        if (!IS_SOLID(fields->Label(ci + 1, cj + 1))) {
            OMP_PRAGMA(omp atomic)
            fields->countAliveParticles->A[fields->countAliveParticles->nx * cj + ci]
                += varType(1);
        }
    }
}
