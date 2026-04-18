#include "PIC.hpp"
#include <random>

varType PIC::rand01() {
  static std::mt19937 rng(std::random_device{}());
  static std::uniform_real_distribution<varType> dist(varType(0), varType(1));
  return dist(rng);
}

// specific interpolation to U-type grid
varType PIC::interpolateU(const Grid2D &g, const varType x,
                          const varType y) const {
  const varType i_real = x / dx;
  const varType j_real = y / dy - REAL_LITERAL(0.5);
  int i = static_cast<int>(std::floor(i_real));
  int j = static_cast<int>(std::floor(j_real));
  const varType fx = i_real - varType(i);
  const varType fy = j_real - varType(j);
  i = std::clamp(i, 0, g.nx - 2);
  j = std::clamp(j, 0, g.ny - 2);
  return (1 - fy) * ((1 - fx) * g.Get(i, j) + fx * g.Get(i + 1, j)) +
         fy * ((1 - fx) * g.Get(i, j + 1) + fx * g.Get(i + 1, j + 1));
}

// specific interpolation to V-type grid
varType PIC::interpolateV(const Grid2D &g, const varType x,
                          const varType y) const {
  const varType i_real = x / dx - REAL_LITERAL(0.5);
  const varType j_real = y / dy;
  int i = static_cast<int>(std::floor(i_real));
  int j = static_cast<int>(std::floor(j_real));
  const varType fx = i_real - varType(i);
  const varType fy = j_real - varType(j);
  i = std::clamp(i, 0, g.nx - 2);
  j = std::clamp(j, 0, g.ny - 2);
  return (1 - fy) * ((1 - fx) * g.Get(i, j) + fx * g.Get(i + 1, j)) +
         fy * ((1 - fx) * g.Get(i, j + 1) + fx * g.Get(i + 1, j + 1));
}