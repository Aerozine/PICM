#include "core/Fields.hpp"
#include "core/Grid2D.hpp"
#include "core/OutputWriter.hpp"
#include "core/Parameters.hpp"
#include "core/project.hpp"
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

int main(int argc, char *argv[]) {
  /*
   * Macro allow print debug in some build cases
   */
#ifndef NDEBUG
  // maybe stream selector or sth  ?
  std::cout << "Compiled with debug mode" << std::endl;
#endif
  // man page ?
  Parameters params;
  if (!params.parseCommandLine(argc, argv)) {
    exit(1);
  }
  // using operator overload
  std::cout << params << std::endl;

  Grid2D grid(params.nx, params.ny);

  Fields2D fields(params.nx, params.ny, params.density, params.dt, params.dx, params.dy);
  // fields.u.FillRandom();
  fields.u.InitRectangle(10.0);
  // fields.InitPotentialGradient(1.0, 1, 1);
  // Grid2D uNorm = fields.VelocityNormCenterGrid();
  Project project(fields);

  fields.Div();

  // generate in the folder result and the simulation.pvd file
  OutputWriter uWriter("results", "u");
  OutputWriter pWriter("results", "p");
  OutputWriter vWriter("results", "v");
  OutputWriter divWriter("results", "div");
  // OutputWriter uNormWriter("results", "uNorm");

  // do 10 step to check if everythings works
  const int num_steps = 100;

  // generating random grid data noise
  for (int t = 0; t < num_steps; ++t) {
    for (int iy = 0; iy < params.ny; ++iy) {
      for (int ix = 0; ix < params.nx; ++ix) {
        /*double x = ix * dx;
        double y = iy * dy;

        varType val = std::sin(2.0 * M_PI * (x - 0.1 * t))
                                      * std::cos(2.0 * M_PI * y);
        grid.SET(ix, iy, val);*/

         project.MakeIncompressible();
      }
    }
    if (!uWriter.writeGrid2D(fields.u, "u") ||
        !vWriter.writeGrid2D(fields.v, "v") ||
        !pWriter.writeGrid2D(fields.p, "p") ||
        !divWriter.writeGrid2D(fields.div, "div")) {
      std::cerr << "Failed to write step " << t << std::endl;
      return 1;
    }
  }
  return 0;
}
