#pragma once
#include "Grid2D.hpp"
#include <cstdint>
#include <string>
#include <vector>
#define FIELD_USOLID 0.0

#define IS_BC(cell) (cell & ((1<<1) | (1<<2) | (1<<3) | (1<<4)))
#define IS_FLUID(cell) (!(cell & (1<<0)))
#define IS_SOLID(cell) (cell &(1<<0))
class Fields2D {
public:
  /// @brief Possible states for a grid cell.

  enum CellType : uint16_t {
    FLUID = 0,     
    SOLID = 1 << 0, 
    AIR = 1 << 1,

    BC_U = 1 << 2,
    BC_V = 1 << 3,
    BC_P = 1 << 4,
    BC_S = 1 << 5,

    IC_U = 1 << 6,
    IC_V = 1 << 7,
    IC_P = 1 << 8,
    IC_S = 1 << 9 
  };

  int nx,ny;
  varType density; ///< Fluid density (kg/m³).
  varType dt,dx,dy;
  Grid2D u,v,p,div;
  Grid2D normVelocity;
  Grid2D smokeMap;


    Grid2D u_sum, u_weight;
    Grid2D v_sum, v_weight;

    Grid2D countAliveParticles;
  // TODO improve string method to sth else and maybe particularize ?
  Fields2D(int nx, int ny, varType density, varType dt, varType dx, varType dy,
           std::string method, std::string freeSurface)
      : nx(nx), ny(ny), density(density), dt(dt), dx(dx), dy(dy), u(nx + 1, ny),
        v(nx, ny + 1), p(nx + 2, ny + 2), div(nx, ny), normVelocity(nx, ny),
        smokeMap(nx, ny),
        u_sum(method == "PIC" ? nx + 1 : 0, method == "PIC" ? ny : 0),
        u_weight(method == "PIC" ? nx + 1 : 0, method == "PIC" ? ny : 0),
        v_sum(method == "PIC" ? nx : 0, method == "PIC" ? ny + 1 : 0),
        v_weight(method == "PIC" ? nx : 0, method == "PIC" ? ny + 1 : 0),
        countAliveParticles(method == "PIC" ? nx : 0, method == "PIC" ? ny : 0),
        labels(static_cast<std::size_t>(nx + 2) * (ny + 2), 
        freeSurface == "yes" ? Fields2D::AIR : Fields2D::FLUID) {}

  /// @brief Return the cell type of cell (i, j).
  [[nodiscard]] inline CellType Label(int i, int j) const noexcept {
  return static_cast<CellType>(labels[idx(i, j)]);
  }

  /// @brief Set the cell type of cell (i, j).
  void SetLabel(int i, int j, CellType t) {
    labels[idx(i, j)] |= static_cast<uint16_t>(t);
  }

  ///@brief Compute the divergence to div
  void Div();

  /// @brief interpolate the velocity magnitude per cell to @c normVelocity
  void VelocityNormCenterGrid();

private:
  std::vector<uint16_t> labels; ///< Flat cell-type array, same layout as p.
  // TODO improve string method to sth else and maybe particularize ?

  [[nodiscard]] inline int idx(int i, int j) const noexcept { return nx * j + i; }
};
