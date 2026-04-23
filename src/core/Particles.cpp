#include "Particles.hpp"
#include <random>

// ── thread-safe rand01 ───────────────────────────────────────────────────────
// Each thread gets its own RNG seeded from a shared device, so RefillParticles
// can safely be parallelised in the future without a data race.
varType rand01() {
    // thread_local: one instance per thread, constructed lazily.
    thread_local std::mt19937 rng(std::random_device{}());
    thread_local std::uniform_real_distribution<varType> dist(
        varType(0), varType(1));
    return dist(rng);
}

// ── InitParticleGrid ─────────────────────────────────────────────────────────
void Particles::InitParticleGrid(const Fields2D &fields) {
    pos_x.clear(); pos_y.clear();
    vel_x.clear(); vel_y.clear();
    if (needsAffine) {
        cu_x.clear(); cu_y.clear();
        cv_x.clear(); cv_y.clear();
    }

    // Count fluid cells first so we can reserve exactly once.
    int fluidCells = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            if (IS_FLUID(fields.Label(i + 1, j + 1)))
                ++fluidCells;

    const std::size_t cap =
        static_cast<std::size_t>(fluidCells) * ppcx * ppcy;

    pos_x.reserve(cap); pos_y.reserve(cap);
    vel_x.reserve(cap); vel_y.reserve(cap);
    if (needsAffine) {
        cu_x.reserve(cap); cu_y.reserve(cap);
        cv_x.reserve(cap); cv_y.reserve(cap);
    }

    for (int jcell = 0; jcell < ny; ++jcell) {
        for (int icell = 0; icell < nx; ++icell) {
            if (!IS_FLUID(fields.Label(icell + 1, jcell + 1)))
                continue;

            for (int a = 0; a < ppcx; ++a) {
                for (int b = 0; b < ppcy; ++b) {
                    const varType x = (icell + rand01()) * dx;
                    const varType y = (jcell + rand01()) * dy;

                    pos_x.push_back(x); pos_y.push_back(y);
                    vel_x.push_back(varType(0));
                    vel_y.push_back(varType(0));
                    if (needsAffine) {
                        cu_x.push_back(varType(0)); cu_y.push_back(varType(0));
                        cv_x.push_back(varType(0)); cv_y.push_back(varType(0));
                    }
                }
            }
        }
    }
}
