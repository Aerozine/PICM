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

  // Create OutputWriter
  OutputWriter writer(params.folder, params.filename);
  // Grid2D grid(nx, ny);
  // Grid2D grid = Grid2D::InitRandomGrid(nx,ny);
  
  varType density = 1000;
  Fields2D fields(params.nx, params.ny, density, params.dt, params.dx, params.dy);
  // fields.u.FillRandom();
  fields.u.InitRectangle(-1.0);
  // fields.InitPotentialGradient(1.0, 1, 1);

  Project project(fields);

  fields.Div();

  OutputWriter uWriter("results", "u");
  OutputWriter pWriter("results", "p");
  OutputWriter vWriter("results", "v");
  OutputWriter uNormWriter("results", "uNorm");
  OutputWriter divWriter("results", "div");

  const int num_steps = 20;

  for (int t = 0; t < num_steps; ++t) {
    Grid2D uNorm = fields.VelocityNormCenterGrid(); // faire mieux 
    if (!uWriter.writeGrid2D(fields.u, "u") ||
        !vWriter.writeGrid2D(fields.v, "v") ||
        !uNormWriter.writeGrid2D(uNorm, "uNorm") ||
        !pWriter.writeGrid2D(fields.p, "p") ||
        !divWriter.writeGrid2D(fields.div, "div")) {
      std::cerr << "Failed to write step " << t << std::endl;
      return 1;
    }
    project.MakeIncompressible();

    // Write with actual time value
  }
  return 0;
}
