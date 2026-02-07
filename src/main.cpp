#include "core/Grid2D.hpp"
#include "core/OutputWriter.hpp"
#include "core/Parameters.hpp"
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
typedef float varType;

int main(int argc, char *argv[]) {
  // man page ?
  Parameters params;
  if (!params.parseCommandLine(argc, argv)) {
    return 1;
  }
  // first test
  // std::cout << "Hello, World!" << std::endl;
  // second test
  // Grid2D grid(5,5);
  // std::cout << grid.A << std::endl;
  // grid.A(2,3)=1.0;
  // print matrix
  // std::cout << grid.A << std::endl;

  // params.print();
  // using operator overload
  std::cout << params << std::endl;

  // example of grid to vtk
  //
  const size_t nx = 50;
  const size_t ny = 40;
  // Grid2D grid(nx, ny);
  Grid2D grid = Grid2D::InitRandomGrid(nx,ny);


  // generate in the folder result and the simulation.pvd file
  OutputWriter writer("results", "simulation");
  // do 10 step to check if everythings works
  const int num_steps = 100;
  const double dx = 1.0 / (nx - 1);
  const double dy = 1.0 / (ny - 1);
  // generating random grid data noise
  for (int t = 0; t < num_steps; ++t) {
    for (size_t iy = 0; iy < ny; ++iy) {
      for (size_t ix = 0; ix < nx; ++ix) {
        double x = ix * dx;
        double y = iy * dy;
        
        // varType val = std::sin(2.0 * M_PI * (x - 0.1 * t)) * std::cos(2.0 * M_PI * y);
        // grid.SET(ix, iy, val);
      }
    }
    // write the grid in the
    if (!writer.writeGrid2D(grid, "random")) {
      std::cerr << "Failed to write step " << t << std::endl;
      return 1;
    }
  }

  // Finalise the .pvd
  // destructor would do it too, but being explicit
  writer.finalisePVD();
  return 0;
}
