
/*
#include "APIC.hpp"

APIC::APIC(Parameters &params) : PIC(params) {}

// ── per-cell affine scatter (no atomics, writes into caller-owned buffer) ────
void APIC::scatterToGrid_APIC(varType xg, varType yg,
                               varType xp, varType yp,
                               varType baseVal, varType cX, varType cY,
                               varType faceOffsetX, varType faceOffsetY,
                               varType* sum, varType* weight,
                               int imax, int jmax) {
    const int i0 = static_cast<int>(std::floor(xg));
    const int j0 = static_cast<int>(std::floor(yg));
    const int radius = params.kernelOrder;

    for (int dj = -radius; dj <= radius; ++dj) {
        for (int di = -radius; di <= radius; ++di) {
            const int i = i0 + di;
            const int j = j0 + dj;
            if (i < 0 || i >= imax || j < 0 || j >= jmax) continue;

            const varType w = hat(xg - varType(i)) * hat(yg - varType(j));
            if (w <= varType(0)) continue;

            const varType xFace = (varType(i) + faceOffsetX) * dx;
            const varType yFace = (varType(j) + faceOffsetY) * dy;
            const varType ox    = xFace - xp;
            const varType oy    = yFace - yp;

            const int flat = imax * j + i;
            sum[flat]    += w * (baseVal + cX * ox + cY * oy);
            weight[flat] += w;
        }
    }
}

// ── scatter one particle's u and v (with affine terms) into private buffers ──
void APIC::scatterParticle_APIC(int idx,
                                 varType* us, varType* uw,
                                 varType* vs, varType* vw) {
    const varType x   = particles->GetX(idx);
    const varType y   = particles->GetY(idx);
    const varType up  = particles->GetU(idx);
    const varType vp  = particles->GetV(idx);
    const varType cuX = particles->GetCuX(idx);
    const varType cuY = particles->GetCuY(idx);
    const varType cvX = particles->GetCvX(idx);
    const varType cvY = particles->GetCvY(idx);

    scatterToGrid_APIC(x / dx,             y / dy - varType(0.5),
                       x, y, up, cuX, cuY,
                       varType(0), varType(0.5),
                       us, uw, nx + 1, ny);

    scatterToGrid_APIC(x / dx - varType(0.5), y / dy,
                       x, y, vp, cvX, cvY,
                       varType(0.5), varType(0),
                       vs, vw, nx, ny + 1);
}

// ── full P2G with thread-local buffers (same pattern as PIC) ─────────────────
void APIC::ProjectParticlesOnGrid() {
    const int nu = (nx + 1) * ny;
    const int nv =  nx      * (ny + 1);

    int nthreads = 1;
    OMP_PRAGMA(omp parallel) { OMP_PRAGMA(omp single) {
#ifdef USE_OPENMP
        nthreads = omp_get_num_threads();
#endif
    }}

    std::vector<varType> u_sum_all(nu * nthreads, varType(0));
    std::vector<varType> u_wgt_all(nu * nthreads, varType(0));
    std::vector<varType> v_sum_all(nv * nthreads, varType(0));
    std::vector<varType> v_wgt_all(nv * nthreads, varType(0));

    OMP_PRAGMA(omp parallel for schedule(static))
    for (int idx = 0; idx < particles->size(); ++idx) {
        int tid = 0;
#ifdef USE_OPENMP
        tid = omp_get_thread_num();
#endif
        scatterParticle_APIC(idx,
            u_sum_all.data() + tid * nu, u_wgt_all.data() + tid * nu,
            v_sum_all.data() + tid * nv, v_wgt_all.data() + tid * nv);
    }

    // reduce
    OMP_PRAGMA(omp parallel for schedule(static))
    for (int k = 0; k < nu; ++k) {
        varType s = varType(0), w = varType(0);
        for (int t = 0; t < nthreads; ++t) {
            s += u_sum_all[t * nu + k];
            w += u_wgt_all[t * nu + k];
        }
        fields->u_sum->A[k] = s;
        fields->u_weight->A[k] = w;
    }
    OMP_PRAGMA(omp parallel for schedule(static))
    for (int k = 0; k < nv; ++k) {
        varType s = varType(0), w = varType(0);
        for (int t = 0; t < nthreads; ++t) {
            s += v_sum_all[t * nv + k];
            w += v_wgt_all[t * nv + k];
        }
        fields->v_sum->A[k] = s;
        fields->v_weight->A[k] = w;
    }

    // normalise — identical to PIC, copied here since ProjectParticlesOnGrid
    // is overridden and PIC's version is not called
    OMP_PRAGMA(omp parallel for collapse(2))
    for (int j = 0; j < fields->u.ny; ++j) {
        for (int i = 0; i < fields->u.nx; ++i) {
            const labeltype left  = fields->Label(i,     j + 1);
            const labeltype right = fields->Label(i + 1, j + 1);
            if (IS_SOLID(left) || IS_SOLID(right)) { fields->u.Set(i, j, varType(0)); continue; }
            if (IS_BC_U(left)) continue;
            const varType w = fields->u_weight->Get(i, j);
            fields->u.Set(i, j, w >= REAL_EPSILON ? fields->u_sum->Get(i, j) / w : varType(0));
        }
    }
    OMP_PRAGMA(omp parallel for collapse(2))
    for (int j = 0; j < fields->v.ny; ++j) {
        for (int i = 0; i < fields->v.nx; ++i) {
            const labeltype bot = fields->Label(i + 1, j);
            const labeltype top = fields->Label(i + 1, j + 1);
            if (IS_SOLID(bot) || IS_SOLID(top)) { fields->v.Set(i, j, varType(0)); continue; }
            if (IS_BC_V(bot)) continue;
            const varType w = fields->v_weight->Get(i, j);
            fields->v.Set(i, j, w >= REAL_EPSILON ? fields->v_sum->Get(i, j) / w : varType(0));
        }
    }
}

// ── G2P unchanged from your original ─────────────────────────────────────────
void APIC::ProjectGridOnParticles() {
    // ... your existing implementation unchanged
}

*/