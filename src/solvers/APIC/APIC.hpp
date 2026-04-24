/*
#pragma once
#include "../PIC/PIC.hpp"

class APIC : public PIC {
public:
    explicit APIC(Parameters &params);

protected:
    // step-level overrides — fine as virtual, called once per step
    void ProjectParticlesOnGrid() override;
    void ProjectGridOnParticles() override;

private:
    // private helpers — never virtual, called per-particle inside this class only
    void scatterToGrid_APIC(varType xg, varType yg,
                             varType xp, varType yp,
                             varType baseVal, varType cX, varType cY,
                             varType faceOffsetX, varType faceOffsetY,
                             varType* sum, varType* weight,
                             int imax, int jmax);

    void scatterParticle_APIC(int idx,
                               varType* us, varType* uw,
                               varType* vs, varType* vw);
};

*/