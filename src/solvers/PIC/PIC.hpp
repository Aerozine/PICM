#pragma once
#include "../../core/Particles.hpp"
#include "../../core/Cloud2D.hpp"
#include "../Solver.hpp"
#include <algorithm>
#include <cmath>
#include <memory>

class PIC : public Solver {
public:
    explicit PIC(Parameters &params);

    void Run()  override;
    void Step() override;

protected:
    // Each grid cell owns its own particle container.
    std::unique_ptr<Cloud2D> cloud;

    inline varType hat(varType r) const {
#ifdef USE_SPEED
        return std::max(varType(0), varType(1) - std::abs(r));
#else
        if (varType(-1.5) <= r && r < varType(-0.5))
            return varType(0.5) * (r + varType(1.5)) * (r + varType(1.5));
        if (varType(-0.5) <= r && r < varType(0.5))
            return varType(0.75) - r * r;
        if (varType(0.5)  <= r && r < varType(1.5))
            return varType(0.5) * (varType(1.5) - r) * (varType(1.5) - r);
        return varType(0);
#endif
    }

    inline varType dhat(varType r) const {
#ifdef USE_SPEED
        return static_cast<varType>((r > varType(-1)) && (r < varType(0))) -
               static_cast<varType>((r > varType(0)) && (r < varType(1)));
#else
        if (varType(-1.5) <= r && r < varType(-0.5)) return r + static_cast<varType>(1.5);
        if (varType(-0.5) <= r && r < varType(0.5))  return -static_cast<varType>(2.0) * r;
        if (varType(0.5)  <= r && r < varType(1.5))  return r - static_cast<varType>(1.5);
        return static_cast<varType>(0);
#endif
    }
    void scatterCell(const Particles &cell,
                     Grid2D &u_sum, Grid2D &u_weight,
                     Grid2D &v_sum, Grid2D &v_weight);

    virtual void ProjectParticlesOnGrid();
    virtual void ProjectGridOnParticles();

    void Advect() override;

    void RefillParticles();
    void UpdateCellState() const;

    void WriteOutput(int step) const;

private:
    std::unique_ptr<OutputWriter> particlesWriter;
};
