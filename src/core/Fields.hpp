#pragma once
#include "Grid2D.hpp"
#include "SolverConfig.hpp"
#include <algorithm>
#include <cmath> // NAN
#include <cstdint>
#include <memory>

#define FIELD_USOLID NAN

#define IS_SOLID(cell) ((((cell)) & Fields2D::SOLID) != 0)
#define IS_AIR(cell) ((((cell)) & Fields2D::AIR) != 0)
#define IS_FLUID(cell) (!IS_SOLID(cell) && !IS_AIR(cell))
#define IS_BC(cell)                                                            \
  ((((cell)) &                                                                 \
    (Fields2D::BC_U | Fields2D::BC_V | Fields2D::BC_P | Fields2D::BC_S)) != 0)
#define IS_BC_U(cell) ((((cell)) & Fields2D::BC_U) != 0)
#define IS_BC_V(cell) ((((cell)) & Fields2D::BC_V) != 0)
#define IS_BC_S(cell) ((((cell)) & Fields2D::BC_S) != 0)
#define IS_BC_P(cell) ((((cell)) & Fields2D::BC_P) != 0)

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

  int nx, ny;
  varType density;
  varType dt, dx, dy;

  // Value members — managed by Grid2D's own Rule-of-Five
  Grid2D u, v, p, div, normVelocity;

  // Optional surface-tension grids — null when surfaceTension == false
  std::unique_ptr<Grid2D> phi;
  std::unique_ptr<Grid2D> kappa;
  std::unique_ptr<Grid2D> normalX;
  std::unique_ptr<Grid2D> normalY;
  std::unique_ptr<Grid2D> interface_u;
  std::unique_ptr<Grid2D> interface_v;

  // Cell labels — (nx+2)*(ny+2), same layout as p
  std::unique_ptr<uint16_t[]> Labels;

  Fields2D(int nx_, int ny_, varType density_, varType dt_, varType dx_,
           varType dy_, bool freeSurface = false, bool surfaceTension = false)
      : nx(nx_), ny(ny_), density(density_), dt(dt_), dx(dx_), dy(dy_),
        u(nx_ + 1, ny_), v(nx_, ny_ + 1), p(nx_ + 2, ny_ + 2), div(nx_, ny_),
        normVelocity(nx_, ny_),
        Labels(std::make_unique<uint16_t[]>(
            static_cast<std::size_t>(nx_ + 2) * (ny_ + 2))) {
    const std::size_t labelCount =
        static_cast<std::size_t>(nx_ + 2) * (ny_ + 2);
    const uint16_t initLabel =
        freeSurface ? static_cast<uint16_t>(AIR) : static_cast<uint16_t>(FLUID);
    std::fill(Labels.get(), Labels.get() + labelCount, initLabel);

    if (surfaceTension) {
      phi        = std::make_unique<Grid2D>(nx_, ny_);
      kappa      = std::make_unique<Grid2D>(nx_, ny_);
      normalX    = std::make_unique<Grid2D>(nx_, ny_);
      normalY    = std::make_unique<Grid2D>(nx_, ny_);
      interface_u = std::make_unique<Grid2D>(nx_ + 1, ny_);
      interface_v = std::make_unique<Grid2D>(nx_, ny_ + 1);
    }
  }

  // unique_ptrs and Grid2D members clean up automatically
  ~Fields2D() = default;

  // Non-copyable; move could be added if ever needed
  Fields2D(const Fields2D &)            = delete;
  Fields2D &operator=(const Fields2D &) = delete;

  [[nodiscard]] inline labeltype Label(int i, int j) const noexcept {
    return Labels[idx(i, j)];
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
    Labels[idx(i, j)] = Labels[idx(i, j)] & ~CELL_TYPE_MASK;
  }

  [[nodiscard]] inline int idx(int i, int j) const noexcept {
    return (nx + 2) * j + i;
  }

  void UpdateDivNorm();
  void Div();
  void VelocityNormCenterGrid();
};
