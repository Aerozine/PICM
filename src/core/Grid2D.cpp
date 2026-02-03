#include "Grid2D.hpp"
// initialisation of the structure
Grid2D::Grid2D(size_t nx, size_t ny) : nx(nx), ny(ny) {
  A = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(ny, nx);
}
