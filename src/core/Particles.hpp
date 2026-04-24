#pragma once
#include "Fields.hpp"
#include "Parameters.hpp"
#include "Precision.hpp"
#include <cstring>
#include <vector>

class Particles {
public:
  int nx, ny, ppcx, ppcy;
  varType dx, dy;

  /// True when cu/cv are actually stored (APIC, mixed FLIP/affine, …).
  const bool needsAffine;

  std::vector<varType> pos_x, pos_y;
  std::vector<varType> vel_x, vel_y;
  std::vector<varType> cu_x, cu_y; ///< empty when !needsAffine
  std::vector<varType> cv_x, cv_y; ///< empty when !needsAffine

  explicit Particles(Parameters &params)
      : nx(params.nx), ny(params.ny), ppcx(params.ppcx), ppcy(params.ppcy),
        dx(params.dx), dy(params.dy),
        needsAffine(params.solver.method != SolverConfig::Method::VanillaPIC &&
                    params.solver.method != SolverConfig::Method::FLIP &&
                    params.solver.method != SolverConfig::Method::SL) {}

  [[nodiscard]] int size() const noexcept {
    return static_cast<int>(pos_x.size());
  }

  /// Add a particle.  cuX/cuY/cvX/cvY are ignored when !needsAffine.
  void Add(varType x, varType y, varType u, varType v,
           [[maybe_unused]] unsigned /*id*/, varType cuX = varType(0),
           varType cuY = varType(0), varType cvX = varType(0),
           varType cvY = varType(0)) {
    pos_x.push_back(x);
    pos_y.push_back(y);
    vel_x.push_back(u);
    vel_y.push_back(v);
    if (needsAffine) {
      cu_x.push_back(cuX);
      cu_y.push_back(cuY);
      cv_x.push_back(cvX);
      cv_y.push_back(cvY);
    }
  }

  void Remove(int i) {
    const int last = size() - 1;
    if (i != last) {
      pos_x[i] = pos_x[last];
      pos_y[i] = pos_y[last];
      vel_x[i] = vel_x[last];
      vel_y[i] = vel_y[last];
      if (needsAffine) {
        cu_x[i] = cu_x[last];
        cu_y[i] = cu_y[last];
        cv_x[i] = cv_x[last];
        cv_y[i] = cv_y[last];
      }
    }
    pos_x.pop_back();
    pos_y.pop_back();
    vel_x.pop_back();
    vel_y.pop_back();
    if (needsAffine) {
      cu_x.pop_back();
      cu_y.pop_back();
      cv_x.pop_back();
      cv_y.pop_back();
    }
  }

  int Compact(const std::vector<uint8_t> &keep) {
    const int n = size();
    int write = 0;
    for (int read = 0; read < n; ++read) {
      if (!keep[read])
        continue;
      if (read != write) {
        pos_x[write] = pos_x[read];
        pos_y[write] = pos_y[read];
        vel_x[write] = vel_x[read];
        vel_y[write] = vel_y[read];
        if (needsAffine) {
          cu_x[write] = cu_x[read];
          cu_y[write] = cu_y[read];
          cv_x[write] = cv_x[read];
          cv_y[write] = cv_y[read];
        }
      }
      ++write;
    }
    resize(write);
    return write;
  }

  [[nodiscard]] varType GetX(int i) const noexcept { return pos_x[i]; }
  [[nodiscard]] varType GetY(int i) const noexcept { return pos_y[i]; }
  [[nodiscard]] varType GetU(int i) const noexcept { return vel_x[i]; }
  [[nodiscard]] varType GetV(int i) const noexcept { return vel_y[i]; }
  [[nodiscard]] varType GetCuX(int i) const noexcept {
    return needsAffine ? cu_x[i] : varType(0);
  }
  [[nodiscard]] varType GetCuY(int i) const noexcept {
    return needsAffine ? cu_y[i] : varType(0);
  }
  [[nodiscard]] varType GetCvX(int i) const noexcept {
    return needsAffine ? cv_x[i] : varType(0);
  }
  [[nodiscard]] varType GetCvY(int i) const noexcept {
    return needsAffine ? cv_y[i] : varType(0);
  }

  void SetX(int i, varType v) noexcept { pos_x[i] = v; }
  void SetY(int i, varType v) noexcept { pos_y[i] = v; }
  void SetU(int i, varType v) noexcept { vel_x[i] = v; }
  void SetV(int i, varType v) noexcept { vel_y[i] = v; }
  void SetCu(int i, varType cx, varType cy) noexcept {
    if (needsAffine) {
      cu_x[i] = cx;
      cu_y[i] = cy;
    }
  }
  void SetCv(int i, varType cx, varType cy) noexcept {
    if (needsAffine) {
      cv_x[i] = cx;
      cv_y[i] = cy;
    }
  }

  void InitParticleGrid(const Fields2D &fields);

  [[nodiscard]] varType *PosX() noexcept { return pos_x.data(); }
  [[nodiscard]] const varType *PosX() const noexcept { return pos_x.data(); }
  [[nodiscard]] varType *PosY() noexcept { return pos_y.data(); }
  [[nodiscard]] const varType *PosY() const noexcept { return pos_y.data(); }
  [[nodiscard]] varType *VelX() noexcept { return vel_x.data(); }
  [[nodiscard]] const varType *VelX() const noexcept { return vel_x.data(); }
  [[nodiscard]] varType *VelY() noexcept { return vel_y.data(); }
  [[nodiscard]] const varType *VelY() const noexcept { return vel_y.data(); }

private:
  void resize(int n) {
    pos_x.resize(n);
    pos_y.resize(n);
    vel_x.resize(n);
    vel_y.resize(n);
    if (needsAffine) {
      cu_x.resize(n);
      cu_y.resize(n);
      cv_x.resize(n);
      cv_y.resize(n);
    }
  }
};

varType rand01();
