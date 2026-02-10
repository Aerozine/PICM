#include "Grid2D.hpp"
#include <cassert>
#include <random>

// prefer some macros ?


// utility functions
bool Grid2D::InBounds(size_t i, size_t j) { return (i < nx) && (j < ny); }

void Grid2D::FillRandom() {
  std::mt19937 gen(42); // fixed seed
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  for (size_t j = 0; j < ny; j++) {
    for (size_t i = 0; i < nx; i++) {
      this->Set(i, j, dist(gen));
    }
  }
  return;
}

/*
Grid2D::Real Grid2D::Interpolate(Real x, Real y, Real dx, Real dy){
  size_t i = std::floor(x / dx);
  size_t j = std::floor(y / dy);

  Real restX = x % i;
  Real restY = y % j;

  return 0.0;
}
*/

Real Grid2D::InterpolateCenter(size_t iCenter, size_t jCenter,
                               const size_t field) {
  // iCenter/jCenter = numbering of pressure grid
  Real interpol;

  switch (field) {
  case 0:
    interpol =
        (this->Get(iCenter, jCenter) + this->Get(iCenter, jCenter + 1)) / 2;
    break;

  case 1:
    interpol =
        (this->Get(iCenter, jCenter) + this->Get(iCenter + 1, jCenter)) / 2;
    break;

  default:
    assert(false);
    break;
  }
  return interpol;
}
