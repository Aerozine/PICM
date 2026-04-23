/*
#include "IterativeMethods.hpp"
#include "Utils.hpp"
#include "../Precision.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
// This code is provided only as example , must be implemented by ourself
//@todo why not using Eigen >_<
//@todo sparsify pressure etc
// @todo make it work
struct SparseA {
    // CSR storage more optimized when free surf
    std::vector<int>    row_ptr;  // size: nrows+1
    std::vector<int>    col_id;   // flat p-space column indices
    std::vector<varType> val;      // corresponding values (+4 or -1)
    // Maps sparse row index → flat p-space id
    std::vector<int>    row_id;
    int                 nrows = 0;
};

static inline int PID(int pnx, int i, int j) { return pnx * j + i; }
static inline bool IS_FLUID_CELL(labeltype l) {
    return !IS_SOLID(l) && !IS_AIR(l) && !IS_BC_P(l);
}

// Connectivity: true when neighbour contributes a -1 off-diagonal entry.
// Matches neighbourSumVec / gsUpdate exactly.
static inline bool fluidLeft  (const Fields2D &f, int i, int j) {
    const labeltype nb = f.Label(i-1, j);
    return !IS_SOLID(nb) && !IS_BC_U(nb) && !IS_AIR(nb);
}
static inline bool fluidRight (const Fields2D &f, int i, int j) {
    // Neumann on right face uses cur label (same asymmetry as gsUpdate)
    return !IS_SOLID(f.Label(i+1,j)) && !IS_BC_U(f.Label(i,j)) && !IS_AIR(f.Label(i+1,j));
}
static inline bool fluidBottom(const Fields2D &f, int i, int j) {
    const labeltype nb = f.Label(i, j-1);
    return !IS_SOLID(nb) && !IS_BC_V(nb) && !IS_AIR(nb);
}
static inline bool fluidTop   (const Fields2D &f, int i, int j) {
    // Neumann on top face uses cur label (same asymmetry as gsUpdate)
    return !IS_SOLID(f.Label(i,j+1)) && !IS_BC_V(f.Label(i,j)) && !IS_AIR(f.Label(i,j+1));
}

static SparseA buildSparseA(const Fields2D &f, int pnx, int pny) {
    SparseA A;

    // First pass: collect fluid cell ids in row-major order.
    for (int j = 1; j < pny - 1; ++j)
        for (int i = 1; i < pnx - 1; ++i)
            if (IS_FLUID_CELL(f.Label(i, j)))
                A.row_id.push_back(PID(pnx, i, j));

    A.nrows = static_cast<int>(A.row_id.size());
    A.row_ptr.resize(A.nrows + 1);

    // Second pass: build CSR.  Diagonal first, then neighbours in id order.
    int nnz = 0;
    for (int row = 0; row < A.nrows; ++row) {
        A.row_ptr[row] = nnz;

        const int id = A.row_id[row];
        const int i  = id % pnx;
        const int j  = id / pnx;

        // Count actual neighbors to set diagonal correctly
        int diag = 0;
        if (fluidBottom(f, i, j)) ++diag;
        if (fluidLeft  (f, i, j)) ++diag;
        if (fluidRight (f, i, j)) ++diag;
        if (fluidTop   (f, i, j)) ++diag;
        
        // Diagonal = number of fluid neighbors (matches neighbourSum)
        const varType diagonal = static_cast<varType>(diag);
        A.col_id.push_back(id);
        A.val   .push_back(diagonal);
        ++nnz;

        // Off-diagonals = -1, stored in ascending column order.
        if (fluidBottom(f, i, j)) { A.col_id.push_back(id - pnx); A.val.push_back(-1.0); ++nnz; }
        if (fluidLeft  (f, i, j)) { A.col_id.push_back(id - 1  ); A.val.push_back(-1.0); ++nnz; }
        if (fluidRight (f, i, j)) { A.col_id.push_back(id + 1  ); A.val.push_back(-1.0); ++nnz; }
        if (fluidTop   (f, i, j)) { A.col_id.push_back(id + pnx); A.val.push_back(-1.0); ++nnz; }
    }
    A.row_ptr[A.nrows] = nnz;

    return A;
}

// y = A·x  (sparse, fluid cells only; y is zero-initialised before call)
static void applyA(const SparseA &A, const std::vector<varType> &x,
                   std::vector<varType> &y) {
    OMP_PRAGMA(omp parallel for schedule(static))
    for (int row = 0; row < A.nrows; ++row) {
        varType acc = 0.0;
        for (int k = A.row_ptr[row]; k < A.row_ptr[row + 1]; ++k)
            acc += A.val[k] * x[A.col_id[k]];
        y[A.row_id[row]] = acc;
    }
}

//  Common RHS / warm-start / write-back helpers
static varType buildRHS(const Fields2D &f, int pnx, int pny, varType coef,
                       const SparseA &A, std::vector<varType> &b) {
    std::fill(b.begin(), b.end(), 0.0);
    varType b_norm = 0.0;
    for (int row = 0; row < A.nrows; ++row) {
        const int id = A.row_id[row];
        const int i  = id % pnx, j = id / pnx;
        const varType val = -coef * static_cast<varType>(f.div.Get(i-1, j-1));
        b[id]  = val;
        b_norm = std::max(b_norm, std::abs(val));
    }
    return b_norm;
}

static void warmStart(const Fields2D &f, int pnx,
                      const SparseA &A, std::vector<varType> &p) {
    for (int row = 0; row < A.nrows; ++row) {
        const int id = A.row_id[row];
        p[id] = static_cast<varType>(f.p.Get(id % pnx, id / pnx));
    }
}

static void writeBack(Fields2D &f, int pnx,
                      const SparseA &A, const std::vector<varType> &p) {
    for (int row = 0; row < A.nrows; ++row) {
        const int id = A.row_id[row];
        f.p.Set(id % pnx, id / pnx, static_cast<varType>(p[id]));
    }
}

//  CG  (unpreconditioned)
bool solveCG(Fields2D &fields, varType coef, varType beta,
             int maxIters, varType tol) {
    fields.Div();

    const int pnx = fields.p.nx;
    const int pny = fields.p.ny;
    const int N   = pnx * pny;

    const SparseA A = buildSparseA(fields, pnx, pny);

    std::vector<varType> b(N, 0.0), p(N, 0.0);
    const varType b_norm = buildRHS(fields, pnx, pny, coef, A, b);

    if (b_norm < 1e-30) {
#ifndef NDEBUG
        std::cout << "  CG: RHS negligible, skipping.\n";
#endif
        return true;
    }

    warmStart(fields, pnx, A, p);

    // r = b - A*p0
    std::vector<varType> Ap(N, 0.0);
    applyA(A, p, Ap);

    std::vector<varType> r(N), d(N);
    varType sigma = 0.0;
    OMP_PRAGMA(omp parallel for reduction(+:sigma))
    for (int k = 0; k < N; ++k) {
        r[k]  = b[k] - Ap[k];
        d[k]  = r[k];
        sigma += r[k] * r[k];
    }

    varType sigma0 = sigma;
    if (sigma0 < 1e-60) {
#ifndef NDEBUG
        std::cout << "  CG: initial residual negligible.\n";
#endif
        writeBack(fields, pnx, A, p);
        return true;
    }

    std::vector<varType> Ad(N, 0.0);
    bool converged = false;

    for (int it = 0; it < maxIters; ++it) {
        applyA(A, d, Ad);

        varType dAd = 0.0;
        OMP_PRAGMA(omp parallel for reduction(+:dAd))
        for (int k = 0; k < N; ++k)
            dAd += d[k] * Ad[k];

        if (std::abs(dAd) < 1e-60) break;
        const varType alpha = sigma / dAd;

        varType r_inf = 0.0;
        OMP_PRAGMA(omp parallel for reduction(max:r_inf))
        for (int k = 0; k < N; ++k) {
            p[k] += alpha * d[k];
            r[k] -= alpha * Ad[k];
            r_inf = std::max(r_inf, std::abs(r[k]));
        }

        if (r_inf / b_norm <= tol) {
#ifndef NDEBUG
            std::cout << "  CG converged in " << it + 1
                      << " iters, rel.res = " << r_inf / b_norm << '\n';
#endif
            converged = true;
            break;
        }

        varType sigma_new = 0.0;
        OMP_PRAGMA(omp parallel for reduction(+:sigma_new))
        for (int k = 0; k < N; ++k)
            sigma_new += r[k] * r[k];

        const varType beta_cg = sigma_new / sigma;
        sigma = sigma_new;

        OMP_PRAGMA(omp parallel for)
        for (int k = 0; k < N; ++k)
            d[k] = r[k] + beta_cg * d[k];
    }

#ifndef NDEBUG
    if (!converged)
        std::cout << "  CG: reached maxIters = " << maxIters << '\n';
#endif

    writeBack(fields, pnx, A, p);
    return converged;
}

//  MICCG0  (Modified Incomplete Cholesky preconditioned CG)
//
//  MIC(0) factorisation follows Bridson "Fluid Simulation for Computer
//  Graphics" §4.4, adapted to the label-based connectivity above.
//  The matrix A has constant diagonal = 4 for fluid cells.
static void buildPrecon(const Fields2D &f, int pnx, int pny,
                        const SparseA &A, std::vector<varType> &precon) {
    constexpr varType tau   = 0.97;   // fill-in strength (Bridson §4.4)
    constexpr varType sigma = 0.25;   // safety factor to prevent negative pivots
    constexpr varType Adiag = 4.0;

    std::fill(precon.begin(), precon.end(), 0.0);

    for (int row = 0; row < A.nrows; ++row) {
        const int id = A.row_id[row];
        const int i  = id % pnx, j = id / pnx;
        varType e = Adiag;

        // Left neighbour (id-1): already processed
        if (fluidLeft(f, i, j)) {
            const varType t = precon[id - 1];
            e -= t * t;
            if (fluidTop(f, i-1, j))  // cross fill-in
                e -= tau * t * t;
        }

        // Bottom neighbour (id-pnx): already processed
        if (fluidBottom(f, i, j)) {
            const varType t = precon[id - pnx];
            e -= t * t;
            if (fluidRight(f, i, j-1))  // cross fill-in
                e -= tau * t * t;
        }

        if (e < sigma * Adiag) e = Adiag;   // safety clamp
        precon[id] = 1.0 / std::sqrt(e);
    }
}

// z = M^-1 r  via forward (L q = r) then backward (L^T z = q) solve.
static void applyPrecon(const Fields2D &f, int pnx, int pny,
                        const SparseA &A,
                        const std::vector<varType> &precon,
                        const std::vector<varType> &r,
                        std::vector<varType>       &z,
                        std::vector<varType>       &q) {
    std::fill(q.begin(), q.end(), 0.0);
    std::fill(z.begin(), z.end(), 0.0);

    // Forward: L q = r  (row-major = lower-triangular order)
    for (int row = 0; row < A.nrows; ++row) {
        const int id = A.row_id[row];
        const int i  = id % pnx, j = id / pnx;
        varType t = r[id];
        if (fluidLeft  (f, i, j)) t -= precon[id - 1  ] * q[id - 1  ];
        if (fluidBottom(f, i, j)) t -= precon[id - pnx] * q[id - pnx];
        q[id] = t * precon[id];
    }

    // Backward: L^T z = q  (reverse row-major = upper-triangular order)
    for (int row = A.nrows - 1; row >= 0; --row) {
        const int id = A.row_id[row];
        const int i  = id % pnx, j = id / pnx;
        varType t = q[id];
        if (fluidRight(f, i, j)) t -= precon[id] * z[id + 1  ];
        if (fluidTop  (f, i, j)) t -= precon[id] * z[id + pnx];
        z[id] = t * precon[id];
    }
}

bool solveMICCG0(Fields2D &fields, varType coef, int maxIters, varType tol) {
    fields.Div();

    const int pnx = fields.p.nx;
    const int pny = fields.p.ny;
    const int N   = pnx * pny;

    const SparseA A = buildSparseA(fields, pnx, pny);

    std::vector<varType> b(N, 0.0), p(N, 0.0);
    const varType b_norm = buildRHS(fields, pnx, pny, coef, A, b);

    if (b_norm < 1e-30) {
#ifndef NDEBUG
        std::cout << "  MICCG0: RHS negligible, skipping.\n";
#endif
        return true;
    }

    warmStart(fields, pnx, A, p);

    // r = b - A*p0
    std::vector<varType> Ap(N, 0.0);
    applyA(A, p, Ap);

    std::vector<varType> r(N, 0.0);
    varType r0_inf = 0.0;
    OMP_PRAGMA(omp parallel for reduction(max:r0_inf))
    for (int k = 0; k < N; ++k) {
        r[k]   = b[k] - Ap[k];
        r0_inf = std::max(r0_inf, std::abs(r[k]));
    }

    if (r0_inf < 1e-30) {
#ifndef NDEBUG
        std::cout << "  MICCG0: initial residual negligible.\n";
#endif
        writeBack(fields, pnx, A, p);
        return true;
    }

    // Build preconditioner
    std::vector<varType> precon(N, 0.0);
    buildPrecon(fields, pnx, pny, A, precon);

#ifndef NDEBUG
    {
        int zeroPre = 0;
        for (int row = 0; row < A.nrows; ++row)
            if (precon[A.row_id[row]] == 0.0) ++zeroPre;
        if (zeroPre > 0)
            std::cout << "  MICCG0 WARNING: " << zeroPre
                      << " fluid cells have zero preconditioner\n";
    }
#endif

    // PCG loop
    std::vector<varType> z(N, 0.0), s(N, 0.0), tmp(N, 0.0), q(N, 0.0);

    applyPrecon(fields, pnx, pny, A, precon, r, z, q);
    s = z;

    varType sigma = 0.0;
    OMP_PRAGMA(omp parallel for reduction(+:sigma))
    for (int k = 0; k < N; ++k)
        sigma += z[k] * r[k];

#ifndef NDEBUG
    if (sigma <= 0.0)
        std::cout << "  MICCG0 WARNING: initial sigma = " << sigma
                  << " (preconditioner may not be SPD)\n";
#endif

    bool converged = false;

    for (int it = 0; it < maxIters; ++it) {
        applyA(A, s, tmp);

        varType dot_s_tmp = 0.0;
        OMP_PRAGMA(omp parallel for reduction(+:dot_s_tmp))
        for (int k = 0; k < N; ++k)
            dot_s_tmp += tmp[k] * s[k];

        if (std::abs(dot_s_tmp) < 1e-60) break;
        const varType alpha = sigma / dot_s_tmp;

        varType r_inf = 0.0;
        OMP_PRAGMA(omp parallel for reduction(max:r_inf))
        for (int k = 0; k < N; ++k) {
            p[k] += alpha * s[k];
            r[k] -= alpha * tmp[k];
            r_inf = std::max(r_inf, std::abs(r[k]));
        }

        if (r_inf / r0_inf <= tol) {
#ifndef NDEBUG
            std::cout << "  MICCG0 converged in " << it + 1
                      << " iters, rel.res = " << r_inf / r0_inf << '\n';
#endif
            converged = true;
            break;
        }

        applyPrecon(fields, pnx, pny, A, precon, r, z, q);

        varType sigma_new = 0.0;
        OMP_PRAGMA(omp parallel for reduction(+:sigma_new))
        for (int k = 0; k < N; ++k)
            sigma_new += z[k] * r[k];

        const varType beta_pcg = sigma_new / sigma;
        sigma = sigma_new;

        OMP_PRAGMA(omp parallel for)
        for (int k = 0; k < N; ++k)
            s[k] = z[k] + beta_pcg * s[k];
    }

#ifndef NDEBUG
    if (!converged)
        std::cout << "  MICCG0: reached maxIters = " << maxIters << '\n';
#endif

    writeBack(fields, pnx, A, p);
    return converged;
}
*/