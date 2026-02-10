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

  // example of grid to vtk
  //
  const size_t nx = 50;
  const size_t ny = 40;
  const float dx = 0.02;
  const float dy = 0.02;
  const float dt = 0.001;
  const float density = 10;

  Grid2D grid(params.nx, params.ny);

  // Create OutputWriter
  OutputWriter writer(params.folder, params.filename);
  // Grid2D grid(nx, ny);
  // Grid2D grid = Grid2D::InitRandomGrid(nx,ny);

  Fields2D fields(nx, ny, density, dt, dx, dy);
  // fields.u.FillRandom();
  fields.u.InitRectangle(10.0);
  // fields.InitPotentialGradient(1.0, 1, 1);
  // Grid2D uNorm = fields.VelocityNormCenterGrid();

  Project project(fields);

  fields.Div();

  OutputWriter uWriter("results", "u");
  OutputWriter pWriter("results", "p");
  OutputWriter vWriter("results", "v");
  OutputWriter divWriter("results", "div");

  const int num_steps = 10;

  for (int t = 0; t < num_steps; ++t) {

    project.MakeIncompressible();
   
    if (!uWriter.writeGrid2D(fields.u, "u") ||
        !vWriter.writeGrid2D(fields.v, "v") ||
        !pWriter.writeGrid2D(fields.p, "p") ||
        !divWriter.writeGrid2D(fields.div, "div")) {
      std::cerr << "Failed to write step " << t << std::endl;
      return 1;
    }
    // Write with actual time value
  }
  return 0;
}
