#pragma once
// simple Eigen lib dense matrix
#include <Eigen/Dense>
// including features of cpp
#include <array>
#include <vector>
class Grid2D {
public:
  Grid2D(size_t nx, size_t ny);
  size_t nx, ny;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A;
};
