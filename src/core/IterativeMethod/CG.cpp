#include "../Precision.hpp"
#include "IterativeMethods.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#ifdef HAVE_EIGEN
#include <Eigen/Sparse>
#endif

namespace {

bool fallbackToRBGS(Fields2D &fields, varType coef, int maxIters, varType tol,
                    varType beta, [[maybe_unused]] const char *solverName) {
#ifndef NDEBUG
  std::cout << "  " << solverName
            << ": Eigen not detected, falling back to RedBlackGS.\n";
#endif
  solveRedBlackGaussSeidel(fields, fields.nx, fields.ny, coef, maxIters, tol,
                           beta);
  return true;
}

#ifdef HAVE_EIGEN

using SparseMatrix = Eigen::SparseMatrix<varType, Eigen::RowMajor>;
using Triplet = Eigen::Triplet<varType>;
using Vector = Eigen::Matrix<varType, Eigen::Dynamic, 1>;

struct PressureSystem {
  int pnx = 0;
  int pny = 0;
  std::vector<int> rowOfCell;
  std::vector<int> cellOfRow;
  SparseMatrix A;
};

struct IC0Preconditioner {
  struct Entry {
    int col;
    varType value;
  };

  std::vector<std::vector<Entry>> lower;
  std::vector<std::vector<Entry>> upper;
  std::vector<varType> diag;
  int clamped = 0;
};

[[nodiscard]] inline int pressureId(const int pnx, const int i,
                                    const int j) noexcept {
  return pnx * j + i;
}

[[nodiscard]] inline bool isUnknownPressureCell(labeltype lbl) noexcept {
  return !IS_SOLID(lbl) && !IS_AIR(lbl) && !IS_BC_P(lbl);
}

[[nodiscard]] inline bool isKnownPressureCell(labeltype lbl) noexcept {
  return !IS_SOLID(lbl) && !IS_AIR(lbl) && IS_BC_P(lbl);
}

[[nodiscard]] bool buildPressureSystem(const Fields2D &fields, varType coef,
                                       PressureSystem &system, Vector &b,
                                       Vector &x0) {
  system.pnx = fields.p.nx;
  system.pny = fields.p.ny;
  system.rowOfCell.assign(static_cast<std::size_t>(system.pnx) * system.pny,
                          -1);
  system.cellOfRow.clear();

  for (int j = 1; j < system.pny - 1; ++j) {
    for (int i = 1; i < system.pnx - 1; ++i) {
      if (!isUnknownPressureCell(fields.Label(i, j)))
        continue;

      const int pid = pressureId(system.pnx, i, j);
      system.rowOfCell[pid] = static_cast<int>(system.cellOfRow.size());
      system.cellOfRow.push_back(pid);
    }
  }

  const int n = static_cast<int>(system.cellOfRow.size());
  system.A.resize(n, n);
  b = Vector::Zero(n);
  x0 = Vector::Zero(n);

  if (n == 0)
    return true;

  std::vector<Triplet> triplets;
  triplets.reserve(static_cast<std::size_t>(n) * 5);

  for (int row = 0; row < n; ++row) {
    const int pid = system.cellOfRow[row];
    const int i = pid % system.pnx;
    const int j = pid / system.pnx;
    const labeltype cur = fields.Label(i, j);

    varType diag = varType(4);
    varType rhs = -coef * fields.div.Get(i - 1, j - 1);

    const auto addUnknownNeighbor = [&](const int ni, const int nj) {
      const int nid = pressureId(system.pnx, ni, nj);
      const int col = system.rowOfCell[nid];
      if (col >= 0)
        triplets.emplace_back(row, col, -varType(1));
    };

    // Left face / neighbour.
    {
      const labeltype nb = fields.Label(i - 1, j);
      if (IS_SOLID(nb) || IS_BC_U(nb)) {
        diag -= varType(1);
      } else if (isKnownPressureCell(nb)) {
        rhs += fields.p.Get(i - 1, j);
      } else if (!IS_AIR(nb)) {
        addUnknownNeighbor(i - 1, j);
      }
    }

    // Right face / neighbour.
    {
      const labeltype nb = fields.Label(i + 1, j);
      if (IS_SOLID(nb) || IS_BC_U(cur)) {
        diag -= varType(1);
      } else if (isKnownPressureCell(nb)) {
        rhs += fields.p.Get(i + 1, j);
      } else if (!IS_AIR(nb)) {
        addUnknownNeighbor(i + 1, j);
      }
    }

    // Bottom face / neighbour.
    {
      const labeltype nb = fields.Label(i, j - 1);
      if (IS_SOLID(nb) || IS_BC_V(nb)) {
        diag -= varType(1);
      } else if (isKnownPressureCell(nb)) {
        rhs += fields.p.Get(i, j - 1);
      } else if (!IS_AIR(nb)) {
        addUnknownNeighbor(i, j - 1);
      }
    }

    // Top face / neighbour.
    {
      const labeltype nb = fields.Label(i, j + 1);
      if (IS_SOLID(nb) || IS_BC_V(cur)) {
        diag -= varType(1);
      } else if (isKnownPressureCell(nb)) {
        rhs += fields.p.Get(i, j + 1);
      } else if (!IS_AIR(nb)) {
        addUnknownNeighbor(i, j + 1);
      }
    }

    if (!std::isfinite(diag) || diag <= REAL_EPSILON)
      return false;

    triplets.emplace_back(row, row, diag);
    b(row) = rhs;
    x0(row) = fields.p.Get(i, j);
  }

  system.A.setFromTriplets(triplets.begin(), triplets.end());
  return true;
}

void writeBackPressure(Fields2D &fields, const PressureSystem &system,
                       const Vector &x) {
  for (int row = 0; row < x.size(); ++row) {
    const int pid = system.cellOfRow[row];
    fields.p.Set(pid % system.pnx, pid / system.pnx, x(row));
  }
}

[[nodiscard]] inline varType safeResidualFloor(varType value) noexcept {
  return std::max(REAL_EPSILON, std::abs(value) * REAL_LITERAL(1e-4));
}

[[nodiscard]] varType
findLowerValue(const std::vector<IC0Preconditioner::Entry> &row,
               const int col) noexcept {
  for (const auto &entry : row) {
    if (entry.col == col)
      return entry.value;
    if (entry.col > col)
      break;
  }
  return varType(0);
}

[[nodiscard]] bool buildIC0(const SparseMatrix &A, IC0Preconditioner &M) {
  const int n = A.rows();
  M.lower.assign(n, {});
  M.upper.assign(n, {});
  M.diag.assign(n, varType(0));
  M.clamped = 0;

  for (int row = 0; row < n; ++row) {
    for (SparseMatrix::InnerIterator it(A, row); it; ++it) {
      if (it.col() < row)
        M.lower[row].push_back({static_cast<int>(it.col()), it.value()});
    }
  }

  for (int row = 0; row < n; ++row) {
    auto &lowerRow = M.lower[row];

    for (std::size_t e = 0; e < lowerRow.size(); ++e) {
      const int col = lowerRow[e].col;
      varType sum = lowerRow[e].value;

      for (std::size_t prev = 0; prev < e; ++prev) {
        const int k = lowerRow[prev].col;
        sum -= lowerRow[prev].value * findLowerValue(M.lower[col], k);
      }

      if (std::abs(M.diag[col]) <= REAL_EPSILON)
        return false;

      lowerRow[e].value = sum / M.diag[col];
    }

    varType d = A.coeff(row, row);
    for (const auto &entry : lowerRow)
      d -= entry.value * entry.value;

    const varType floor = safeResidualFloor(A.coeff(row, row));
    if (!std::isfinite(d))
      return false;
    if (d <= floor) {
      d = floor;
      ++M.clamped;
    }

    M.diag[row] = std::sqrt(d);
  }

  for (int row = 0; row < n; ++row) {
    for (const auto &entry : M.lower[row])
      M.upper[entry.col].push_back({row, entry.value});
  }

  return true;
}

void applyIC0(const IC0Preconditioner &M, const Vector &r, Vector &z) {
  const int n = static_cast<int>(M.diag.size());
  Vector y = Vector::Zero(n);
  z = Vector::Zero(n);

  for (int row = 0; row < n; ++row) {
    varType sum = r(row);
    for (const auto &entry : M.lower[row])
      sum -= entry.value * y(entry.col);
    y(row) = sum / M.diag[row];
  }

  for (int row = n - 1; row >= 0; --row) {
    varType sum = y(row);
    for (const auto &entry : M.upper[row])
      sum -= entry.value * z(entry.col);
    z(row) = sum / M.diag[row];
  }
}

bool solveSparseCGSystem(const PressureSystem &system, const Vector &b,
                         Vector &x, int maxIters, varType tol,
                         [[maybe_unused]] const char *solverName) {
  if (b.size() == 0)
    return true;

  const varType bNorm = b.template lpNorm<Eigen::Infinity>();
  if (bNorm <= REAL_EPSILON)
    return true;

  Vector r = b - system.A * x;
  Vector d = r;
  varType sigma = r.dot(r);

  if (sigma <= REAL_EPSILON)
    return true;

  bool converged = false;
  int it = 0;

  for (; it < maxIters; ++it) {
    const Vector Ad = system.A * d;
    const varType dAd = d.dot(Ad);

    if (!std::isfinite(dAd) || std::abs(dAd) <= REAL_EPSILON)
      break;

    const varType alpha = sigma / dAd;
    x += alpha * d;
    r -= alpha * Ad;

    const varType relRes = r.template lpNorm<Eigen::Infinity>() / bNorm;
    if (relRes <= tol) {
      converged = true;
      ++it;
      debugSolverConverged(solverName, it, relRes);
      break;
    }

    const varType sigmaNew = r.dot(r);
    if (!std::isfinite(sigmaNew) || sigmaNew <= REAL_EPSILON)
      break;

    d = r + (sigmaNew / sigma) * d;
    sigma = sigmaNew;
  }

  if (!converged)
    debugSolverMaxIters(solverName, maxIters);

  return converged;
}

bool solveSparsePCGSystem(const PressureSystem &system,
                          const IC0Preconditioner &precon, const Vector &b,
                          Vector &x, int maxIters, varType tol) {
  if (b.size() == 0)
    return true;

  const varType bNorm = b.template lpNorm<Eigen::Infinity>();
  if (bNorm <= REAL_EPSILON)
    return true;

  Vector r = b - system.A * x;
  Vector z(b.size());
  applyIC0(precon, r, z);

  Vector d = z;
  varType rz = r.dot(z);
  if (!std::isfinite(rz) || std::abs(rz) <= REAL_EPSILON)
    return true;

  bool converged = false;
  int it = 0;

  for (; it < maxIters; ++it) {
    const Vector Ad = system.A * d;
    const varType dAd = d.dot(Ad);

    if (!std::isfinite(dAd) || std::abs(dAd) <= REAL_EPSILON)
      break;

    const varType alpha = rz / dAd;
    x += alpha * d;
    r -= alpha * Ad;

    const varType relRes = r.template lpNorm<Eigen::Infinity>() / bNorm;
    if (relRes <= tol) {
      converged = true;
      ++it;
      debugSolverConverged("MICCG0", it, relRes);
      if (precon.clamped > 0) {
        DBG_PRINTF("MICCG0: IC(0) clamped %d pivot(s)", precon.clamped);
      }
      break;
    }

    applyIC0(precon, r, z);
    const varType rzNew = r.dot(z);
    if (!std::isfinite(rzNew) || std::abs(rzNew) <= REAL_EPSILON)
      break;

    d = z + (rzNew / rz) * d;
    rz = rzNew;
  }

  if (!converged)
    debugSolverMaxIters("MICCG0", maxIters);
  if (precon.clamped > 0) {
    DBG_PRINTF("MICCG0: IC(0) clamped %d pivot(s)", precon.clamped);
  }

  return converged;
}

#endif

} // namespace

bool solveCG(Fields2D &fields, varType coef, varType beta, int maxIters,
             varType tol) {
  fields.Div();

#ifndef HAVE_EIGEN
  return fallbackToRBGS(fields, coef, maxIters, tol, beta, "CG");
#else
  (void)beta;
  PressureSystem system;
  Vector b;
  Vector x;
#ifndef NDEBUG
  const double assemblyStart = GET_TIME();
#endif

  if (!buildPressureSystem(fields, coef, system, b, x)) {
#ifndef NDEBUG
    std::cout << "  CG: sparse matrix assembly produced a non-SPD diagonal, "
                 "falling back to RedBlackGS.\n";
#endif
    return fallbackToRBGS(fields, coef, maxIters, tol, beta, "CG");
  }

#ifndef NDEBUG
  const double assemblyTime = GET_TIME() - assemblyStart;
  const double solveStart = GET_TIME();
#endif
  const bool converged = solveSparseCGSystem(system, b, x, maxIters, tol, "CG");
#ifndef NDEBUG
  DBG_PRINTF("CG: rows=%ld, assembly=%.6fs, solve=%.6fs", system.A.rows(),
             assemblyTime, GET_TIME() - solveStart);
#endif
  writeBackPressure(fields, system, x);
  return converged;
#endif
}

bool solveMICCG0(Fields2D &fields, varType coef, int maxIters, varType tol) {
  fields.Div();

#ifndef HAVE_EIGEN
  return fallbackToRBGS(fields, coef, maxIters, tol, varType(0), "MICCG0");
#else
  PressureSystem system;
  Vector b;
  Vector x;
#ifndef NDEBUG
  const double assemblyStart = GET_TIME();
#endif

  if (!buildPressureSystem(fields, coef, system, b, x)) {
#ifndef NDEBUG
    std::cout << "  MICCG0: sparse matrix assembly produced a non-SPD "
                 "diagonal, falling back to RedBlackGS.\n";
#endif
    return fallbackToRBGS(fields, coef, maxIters, tol, varType(0), "MICCG0");
  }

#ifndef NDEBUG
  const double assemblyTime = GET_TIME() - assemblyStart;
#endif
  IC0Preconditioner precon;
#ifndef NDEBUG
  const double ic0Start = GET_TIME();
#endif
  if (!buildIC0(system.A, precon)) {
#ifndef NDEBUG
    std::cout << "  MICCG0: IC(0) factorisation failed, falling back to CG.\n";
#endif
    const bool converged =
        solveSparseCGSystem(system, b, x, maxIters, tol, "CG");
    writeBackPressure(fields, system, x);
    return converged;
  }

#ifndef NDEBUG
  const double ic0Time = GET_TIME() - ic0Start;
  const double solveStart = GET_TIME();
#endif
  const bool converged =
      solveSparsePCGSystem(system, precon, b, x, maxIters, tol);
#ifndef NDEBUG
  DBG_PRINTF(
      "MICCG0: IC(0) rows=%ld, assembly=%.6fs, factor=%.6fs, solve=%.6fs",
      system.A.rows(), assemblyTime, ic0Time, GET_TIME() - solveStart);
#endif
  writeBackPressure(fields, system, x);
  return converged;
#endif
}
