#pragma once

#ifdef USE_CUDA

#include "../../core/Cloud2D.hpp"
#include "../../core/OutputWriter.hpp"
#include "../../core/Parameters.hpp"
#include "../Solver.hpp"
#include <memory>
#include <vector>

// Full-GPU Particle-In-Cell solver.
// All CUDA device state is hidden behind the PIMPL (DeviceState) so this
// header remains a plain C++ header safe for g++ compilation.
class GPIC : public Solver {
public:
    explicit GPIC(Parameters& params);
    ~GPIC();

    void Run()  override;
    void Step() override;

protected:
    void WriteOutput(int step) const override;

private:
    // PIMPL: defined only in GPIC.cu, compiled by nvcc.
    struct DeviceState;
    std::unique_ptr<DeviceState> dev_;

    // CPU-side particle mirror (kept in sync after every substep).
    // Used for CFL computation, output, and refill.
    mutable std::unique_ptr<Cloud2D> cloud_;

    std::unique_ptr<OutputWriter> particlesWriter_;

    // ── H↔D transfer helpers ─────────────────────────────────────────────
    void UploadParticles();             // cloud_ → device flat SoA
    void DownloadParticles() const;     // device flat SoA → cloud_  (const: modifies mutable cloud_)
    void DownloadAll() const;           // d_u, d_v, d_p, d_labels → fields

    void UploadGrid();                  // fields u, v, Labels → device

    // ── GPU step stages ──────────────────────────────────────────────────
    void GPUProjectParticlesOnGrid();   // P2G
    void GPUMakeIncompressible();       // div → GPU RBGS → pressure gradient
    void GPUProjectGridOnParticles();   // G2P (also applies gravity)
    void GPUAdvect();                   // RK2 + compaction of dead particles
    void GPUCountParticles();           // rebuild per-cell particle counts
    void GPUUpdateCellState();          // update labels + air cleanup (fully GPU)
    void GPURefillParticles();          // refill under-populated cells fully on GPU

    int  computeAdvectionSubstepsGPU() const;  // GPU thrust reduce over particle velocities

    void DownloadLabels() const;         // d_labels → fields
};

#endif // USE_CUDA
