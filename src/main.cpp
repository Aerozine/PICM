#include "core/Fields.hpp"
#include "core/project.hpp"
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
  const float dx = 1.0 / (nx - 1);
  const float dy = 1.0 / (ny - 1);
  const float dt = 0.001;
  const float density = 1000;

  // Grid2D grid(nx, ny);
  // Grid2D grid = Grid2D::InitRandomGrid(nx,ny);
  
  Fields2D fields(nx, ny, density, dt, dx, dy); 
  // fields.u.FillRandom();
  // fields.v.FillRandom();
  fields.InitPotentialGradient(1.0, 1, 1);
  Grid2D uNorm = fields.VelocityNormCenterGrid();

  Project project(fields);
  // project.MakeIncompressible();
  
  fields.Div();
   
  // generate in the folder result and the simulation.pvd file
  OutputWriter uWriter("results", "u");
  OutputWriter vWriter("results", "v");
  OutputWriter divWriter("results", "div");
  OutputWriter uNormWriter("results", "uNorm");

  // do 10 step to check if everythings works
  const int num_steps = 100;
  
  // generating random grid data noise
  for (int t = 0; t < num_steps; ++t) {
    for (size_t iy = 0; iy < ny; ++iy) {
      for (size_t ix = 0; ix < nx; ++ix) {
        /*double x = ix * dx;
        double y = iy * dy;
        
        varType val = std::sin(2.0 * M_PI * (x - 0.1 * t))
                                      * std::cos(2.0 * M_PI * y);
        grid.SET(ix, iy, val);*/
      }
    }
    /*// write the grid in the
    if (!uWriter.writeGrid2D(fields.u, "u") 
      && !vWriter.writeGrid2D(fields.v, "v")
      && !divWriter.writeGrid2D(fields.div, "div")) {
      std::cerr << "Failed to write step " << t << std::endl;
      return 1;
    }*/
    if (!uWriter.writeGrid2D(fields.u, "u") ||
        !vWriter.writeGrid2D(fields.v, "v") ||
        !uNormWriter.writeGrid2D(uNorm, "uNorm") ||
        !divWriter.writeGrid2D(fields.div, "div")) {
      std::cerr << "Failed to write step " << t << std::endl;
      return 1;
    }
  }

  // Finalise the .pvd
  // destructor would do it too, but being explicit
  uWriter.finalisePVD();
  vWriter.finalisePVD();
  uNormWriter.finalisePVD();
  divWriter.finalisePVD();
  return 0;
}
