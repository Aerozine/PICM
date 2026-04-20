#pragma once
#include "Grid2D.hpp"
#include "SolverConfig.hpp"
#include <cstdint>
#include <cstring>  // memset
#include <cstdlib>  // malloc, free, calloc
#include <cmath>    // NAN

#define FIELD_USOLID NAN

#define IS_SOLID(cell)  ((((cell)) & Fields2D::SOLID)  != 0)
#define IS_AIR(cell)    ((((cell)) & Fields2D::AIR)    != 0)
#define IS_FLUID(cell)  (!IS_SOLID(cell) && !IS_AIR(cell))
#define IS_BC(cell)     ((((cell)) & (Fields2D::BC_U | Fields2D::BC_V | \
                                      Fields2D::BC_P | Fields2D::BC_S)) != 0)
#define IS_BC_U(cell)   ((((cell)) & Fields2D::BC_U) != 0)
#define IS_BC_V(cell)   ((((cell)) & Fields2D::BC_V) != 0)
#define IS_BC_S(cell)   ((((cell)) & Fields2D::BC_S) != 0)
#define IS_BC_P(cell)   ((((cell)) & Fields2D::BC_P) != 0)

using labeltype = uint16_t;

class Fields2D {
public:

  enum CellType : uint16_t {
    FLUID = 0,
    SOLID = 1 << 0,
    AIR   = 1 << 1,
    BC_U  = 1 << 2,
    BC_V  = 1 << 3,
    BC_P  = 1 << 4,
    BC_S  = 1 << 5,
    IC_U  = 1 << 6,
    IC_V  = 1 << 7,
    IC_P  = 1 << 8,
    IC_S  = 1 << 9,
  };
  static constexpr uint16_t CELL_TYPE_MASK = SOLID | AIR;

  int     nx, ny;
  varType density;
  varType dt, dx, dy;
  //todo use a template for label ? is there an overhead ?
  Grid2D u, v, p, div, normVelocity, smokeMap;

  Grid2D *u_sum             = nullptr;
  Grid2D *u_weight          = nullptr;
  Grid2D *v_sum             = nullptr;
  Grid2D *v_weight          = nullptr;
  Grid2D *countAliveParticles = nullptr;

  uint16_t *Labels = nullptr;

  // ── constructor ────────────────────────────────────────────────────────────
  // freeSurface: true  → initialise all labels to AIR
  //              false → initialise all labels to FLUID (default, incompressible)
  Fields2D(int nx_, int ny_,
           varType density_, varType dt_, varType dx_, varType dy_,
           const SolverConfig &sol,
           bool freeSurface = false)
      : nx(nx_), ny(ny_),
        density(density_), dt(dt_), dx(dx_), dy(dy_),
        u(nx_ + 1, ny_), v(nx_, ny_ + 1),
        p(nx_ + 2, ny_ + 2), div(nx_, ny_),
        normVelocity(nx_, ny_), smokeMap(nx_, ny_)
  {
    const std::size_t labelCount =
        static_cast<std::size_t>(nx_ + 2) * (ny_ + 2);
    Labels = static_cast<uint16_t *>(
        std::malloc(labelCount * sizeof(uint16_t)));
    const uint16_t initLabel =
        freeSurface ? static_cast<uint16_t>(AIR)
                    : static_cast<uint16_t>(FLUID);
    for (std::size_t k = 0; k < labelCount; ++k)
      Labels[k] = initLabel;

    if (sol.method != SolverConfig::Method::SL) {
      u_sum              = new Grid2D(nx_ + 1, ny_);
      u_weight           = new Grid2D(nx_ + 1, ny_);
      v_sum              = new Grid2D(nx_,     ny_ + 1);
      v_weight           = new Grid2D(nx_,     ny_ + 1);
      countAliveParticles = new Grid2D(nx_,    ny_);
    }
  }

  // ── destructor ─────────────────────────────────────────────────────────────
  ~Fields2D() {
    std::free(Labels);
    delete u_sum;
    delete u_weight;
    delete v_sum;
    delete v_weight;
    delete countAliveParticles;
  }

  // Non-copyable, non-movable (owns raw resources; add if ever needed)
  Fields2D(const Fields2D &) = delete;
  Fields2D &operator=(const Fields2D &) = delete;

  // ── label accessors ────────────────────────────────────────────────────────

  [[nodiscard]] inline CellType Label(int i, int j) const noexcept {
    return static_cast<CellType>(Labels[idx(i, j)]);
  }

  /// OR-in a flag (does not clear existing flags).
  inline void SetLabel(int i, int j, CellType t) {
    Labels[idx(i, j)] |= static_cast<uint16_t>(t);
  }

  /// Replace the label entirely (clears all previous flags).
  inline void ResetLabel(int i, int j, CellType t) {
    Labels[idx(i, j)] = static_cast<uint16_t>(t);
  }

  inline void setSolid(int i, int j) {
    Labels[idx(i, j)] = (Labels[idx(i, j)] & ~CELL_TYPE_MASK) | SOLID;
  }
  inline void setAir(int i, int j) {
    Labels[idx(i, j)] = (Labels[idx(i, j)] & ~CELL_TYPE_MASK) | AIR;
  }
  inline void setFluid(int i, int j) {
    Labels[idx(i, j)] =  Labels[idx(i, j)] & ~CELL_TYPE_MASK;
  }

  [[nodiscard]] inline int idx(int i, int j) const noexcept {
    return (nx + 2) * j + i;
  }

  // ── field computations ─────────────────────────────────────────────────────
  void Div();
  void VelocityNormCenterGrid();
};
