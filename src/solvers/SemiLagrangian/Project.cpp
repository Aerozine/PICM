#include "SemiLagrangian.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

void SemiLagrangian::buildRHS() {
    varType invDx = REAL_LITERAL(1.0) / dx;

    // Initialize RHS to zero
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            rhs[idx(i, j)] = REAL_LITERAL(0.0);
        }
    }

    #pragma omp parallel for collapse(2)
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            if (fields->Label(i, j) != Fields2D::FLUID) {
                continue;
            }

            varType div_val = REAL_LITERAL(0.0);

            // U contribution: (u_right - u_left) / dx
            if (i + 1 < fields->u.nx) {
                div_val += (fields->u.Get(i + 1, j) - fields->u.Get(i, j));
            }

            // V contribution: (v_top - v_bottom) / dy
            if (j + 1 < fields->v.ny) {
                div_val += (fields->v.Get(i, j + 1) - fields->v.Get(i, j));
            }

            rhs[idx(i, j)] = -invDx * div_val;
        }
    }

    #pragma omp parallel for collapse(2)
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            if (fields->Label(i, j) != Fields2D::FLUID)
                continue;

            if (i > 0 && fields->Label(i - 1, j) == Fields2D::SOLID) {
                if (i < fields->u.nx) {
                    rhs[idx(i, j)] += invDx * (fields->u.Get(i, j) - fields->usolid);
                }
            }

            if (i + 1 < nx - 1 && fields->Label(i + 1, j) == Fields2D::SOLID) {
                if (i + 1 < fields->u.nx) {
                    rhs[idx(i, j)] -= invDx * (fields->u.Get(i + 1, j) - fields->usolid);
                }
            }

            if (j > 0 && fields->Label(i, j - 1) == Fields2D::SOLID) {
                if (j < fields->v.ny) {
                    rhs[idx(i, j)] -= invDx * (fields->v.Get(i, j) - fields->usolid);
                }
            }

            if (j + 1 < ny - 1 && fields->Label(i, j + 1) == Fields2D::SOLID) {
                if (j + 1 < fields->v.ny) {
                    rhs[idx(i, j)] += invDx * (fields->v.Get(i, j + 1) - fields->usolid);
                }
            }
        }
    }
}

void SemiLagrangian::buildMatrixA() {
    const varType scaleA = dt / (density * dx * dx);

    // Initialize to zero
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            Adiag[idx(i, j)] = REAL_LITERAL(0.0);
            Ax[idx(i, j)] = REAL_LITERAL(0.0);
            Ay[idx(i, j)] = REAL_LITERAL(0.0);
        }
    }

    // Build matrix coefficients 
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            if (fields->Label(i, j) != Fields2D::FLUID)
                continue;

            varType diag = REAL_LITERAL(0.0);

            // +x neighbor: (i+1, j)
            if (i + 1 < nx - 1 && fields->Label(i + 1, j) != Fields2D::SOLID) {
                diag += scaleA;
                if (fields->Label(i + 1, j) == Fields2D::FLUID) {
                    Ax[idx(i, j)] = -scaleA;
                }
            }

            // -x neighbor: (i-1, j)
            if (i > 0 && fields->Label(i - 1, j) != Fields2D::SOLID) {
                diag += scaleA;
            }

            // +y neighbor: (i, j+1)
            if (j + 1 < ny - 1 && fields->Label(i, j + 1) != Fields2D::SOLID) {
                diag += scaleA;
                if (fields->Label(i, j + 1) == Fields2D::FLUID) {
                    Ay[idx(i, j)] = -scaleA;
                }
            }

            // -y neighbor: (i, j-1)
            if (j > 0 && fields->Label(i, j - 1) != Fields2D::SOLID) {
                diag += scaleA;
            }

            Adiag[idx(i, j)] = diag;
        }
    }
}

varType SemiLagrangian::neighborPressureSum(int i, int j) {
    assert(i < nx - 1 && j < ny - 1);
    varType sum = REAL_LITERAL(0.0);

    if (i + 1 < nx - 1)
        sum += Ax[idx(i, j)] * fields->p.Get(i + 1, j);
    if (i > 0)
        sum += Ax[idx(i - 1, j)] * fields->p.Get(i - 1, j);
    if (j + 1 < ny - 1)
        sum += Ay[idx(i, j)] * fields->p.Get(i, j + 1);
    if (j > 0)
        sum += Ay[idx(i, j - 1)] * fields->p.Get(i, j - 1);

    return sum;
}

// PARALLELIZED JACOBI SOLVER
void SemiLagrangian::solveJacobi(int maxIters, varType tol) {
    Grid2D pNew(nx - 1, ny - 1);
    int iterations = 0;
    for (int it = 0; it < maxIters; it++) {
        varType maxDiff = REAL_LITERAL(0.0);

        // PARALLEL Jacobi iteration
        // reduction max:maxdiff compute maxdiff per thread
        #pragma omp parallel for collapse(2) reduction(max:maxDiff)
        for (int j = 0; j < ny - 1; j++) {
            for (int i = 0; i < nx - 1; i++) {
                if (fields->Label(i, j) != Fields2D::FLUID)
                    continue;

                varType diag = Adiag[idx(i, j)];
                if (diag == REAL_LITERAL(0.0))
                    continue;

                // Compute neighbor sum inline (more efficient than function call)
                varType sum = REAL_LITERAL(0.0);
                if (i + 1 < nx - 1)
                    sum += Ax[idx(i, j)] * fields->p.Get(i + 1, j);
                if (i > 0)
                    sum += Ax[idx(i - 1, j)] * fields->p.Get(i - 1, j);
                if (j + 1 < ny - 1)
                    sum += Ay[idx(i, j)] * fields->p.Get(i, j + 1);
                if (j > 0)
                    sum += Ay[idx(i, j - 1)] * fields->p.Get(i, j - 1);

                varType newVal = (rhs[idx(i, j)] - sum) / diag;

                varType diff = std::abs(newVal - fields->p.Get(i, j));
                if (diff > maxDiff) maxDiff = diff;

                pNew.Set(i, j, newVal);
            }
        }

        // PARALLEL swap p <- pNew
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < ny - 1; j++) {
            for (int i = 0; i < nx - 1; i++) {
                if (fields->Label(i, j) == Fields2D::FLUID) {
                    fields->p.Set(i, j, pNew.Get(i, j));
                }
            }
        }
        iterations = it + 1;
        // Check convergence
        if (maxDiff < tol)
#ifndef NDEBUG
    std::cout << "  Jacobi Max Diff " << tol << "> maxDiff" << std::endl;
#endif
            break;
    }
#ifndef NDEBUG
    std::cout << "  Jacobi converged in " << iterations << " iterations" << std::endl;
#endif
}

void SemiLagrangian::updateVelocities() {
    varType coef = dt / (density * dx);

    // PARALLEL Update U velocities
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < fields->u.ny; ++j) {
        for (int i = 0; i < fields->u.nx; ++i) {
            // Boundary conditions at domain edges
            if (i == 0 || i == fields->u.nx - 1 || j == 0 || j == fields->u.ny - 1) {
                fields->u.Set(i, j, fields->usolid);
                continue;
            }

            // Interior: apply pressure gradient
            if (i > 0 && i - 1 < fields->p.nx && j < fields->p.ny) {
                varType uOld = fields->u.Get(i, j);
                varType pRight = fields->p.Get(i, j);
                varType pLeft = (i > 0) ? fields->p.Get(i - 1, j) : REAL_LITERAL(0.0);
                varType uNew = uOld - coef * (pRight - pLeft);
                fields->u.Set(i, j, uNew);
            }
        }
    }

    // PARALLEL Update V velocities
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < fields->v.ny; ++j) {
        for (int i = 0; i < fields->v.nx; ++i) {
            // Boundary conditions at domain edges
            if (i == 0 || i == fields->v.nx - 1 || j == 0 || j == fields->v.ny - 1) {
                fields->v.Set(i, j, fields->usolid);
                continue;
            }

            // Interior: apply pressure gradient
            if (j > 0 && j - 1 < fields->p.ny && i < fields->p.nx) {
                varType vOld = fields->v.Get(i, j);
                varType pTop = fields->p.Get(i, j);
                varType pBottom = (j > 0) ? fields->p.Get(i, j - 1) : REAL_LITERAL(0.0);
                varType vNew = vOld - coef * (pTop - pBottom);
                fields->v.Set(i, j, vNew);
            }
        }
    }
}

void SemiLagrangian::MakeIncompressible() {
    const int maxIters = 10000;
    const varType tol = REAL_LITERAL(0.0001);

    buildRHS();
    buildMatrixA();
    solveJacobi(maxIters, tol);
    updateVelocities();
}
