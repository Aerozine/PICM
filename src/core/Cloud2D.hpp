#pragma once
#include "Fields.hpp"
#include "Particles.hpp"

#ifdef USE_OPENMP
#include <omp.h>
#endif

// Cloud2D: a nx*ny grid where each cell owns its particles as a small SoA.
// Indexing matches Grid2D row-major storage: cells[nx * j + i] for cell (i, j).
class Cloud2D {
public:
  int nx, ny;
  varType dx, dy;
  std::vector<Particles> cells;

#ifdef USE_OPENMP
  std::vector<omp_lock_t> cellLocks;
#endif

  Cloud2D(Parameters &params)
      : nx(params.nx), ny(params.ny), dx(params.dx), dy(params.dy) {
    cells.reserve(nx * ny);
    for (int k = 0; k < nx * ny; ++k)
      cells.emplace_back(params);

#ifdef USE_OPENMP
    cellLocks.resize(nx * ny);
    for (int k = 0; k < nx * ny; ++k)
      omp_init_lock(&cellLocks[k]);
#endif
  }

  ~Cloud2D() {
#ifdef USE_OPENMP
    for (int k = 0; k < nx * ny; ++k)
      omp_destroy_lock(&cellLocks[k]);
#endif
  }

  Cloud2D(const Cloud2D &) = delete;
  Cloud2D &operator=(const Cloud2D &) = delete;

  [[nodiscard]] inline std::size_t idx(int i, int j) const noexcept {
    return static_cast<std::size_t>(nx) * j + i;
  }

  inline Particles &operator()(int i, int j) { return cells[idx(i, j)]; }
  inline const Particles &operator()(int i, int j) const {
    return cells[idx(i, j)];
  }

  inline int countIn(int i, int j) const { return cells[idx(i, j)].size(); }

  // total particle count across all cells
  int totalSize() const {
    int total = 0;
    for (const auto &c : cells)
      total += c.size();
    return total;
  }

  // seed particles into cells from the fields initial state,
  // mirrors the old Particles::InitParticleGrid
  void InitParticleGrid(const Fields2D &fields, int ppcx, int ppcy) {
    for (auto &c : cells) {
      c.pos_x.clear(), c.pos_y.clear(), c.vel_x.clear(), c.vel_y.clear();
      if (c.needsAffine) {
        c.cu_x.clear(), c.cu_y.clear();
        c.cv_x.clear(), c.cv_y.clear();
      }
    }

    for (int jcell = 0; jcell < ny; ++jcell) {
      for (int icell = 0; icell < nx; ++icell) {
        if (!IS_FLUID(fields.Label(icell + 1, jcell + 1)))
          continue;
        Particles &cell = (*this)(icell, jcell);
        cell.pos_x.reserve(ppcx * ppcy);
        cell.pos_y.reserve(ppcx * ppcy);
        cell.vel_x.reserve(ppcx * ppcy);
        cell.vel_y.reserve(ppcx * ppcy);
        for (int a = 0; a < ppcx; ++a) {
          for (int b = 0; b < ppcy; ++b) {
            const varType x = (icell + rand01()) * dx;
            const varType y = (jcell + rand01()) * dy;
            cell.Add(x, y, varType(0), varType(0), 0);
          }
        }
      }
    }
  }
};
