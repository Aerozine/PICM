#include "core/Grid2D.hpp"
// #include "core/OutputWriter.hpp"
#include "core/BetterOutputWriter.hpp"
#include "core/Parameters.hpp"
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

#ifndef NDEBUG
  std::cout << "Parameters loaded " << params << std::endl;
#endif
  /** if params.solverType == Semilagrangian
   *    Semilagrangian(params)
   *
   *
   */

  Grid2D grid(params.nx, params.ny);

  // Create BetterOutputWriter
  BetterOutputWriter writer(params.folder, params.filename, params.dx,
                            params.dy, 1.0);
  // TODO need to be removed , call and execute solver
  // Time stepping loop
  for (int t = 0; t < 1.0 / params.dt; ++t) {
    // Update grid
    for (int iy = 0; iy < params.ny; ++iy) {
      for (int ix = 0; ix < params.nx; ++ix) {
        double x = ix * params.dx;
        double y = iy * params.dy;
        grid.Set(ix, iy,
                 std::sin(2.0 * M_PI * (x - 0.01 * t)) *
                     std::cos(2.0 * M_PI * y));
      }
    }

    // Write with actual time value
    double time = t * params.dt;
    if (!writer.writeGrid2D(grid, "wave", time)) {
      std::cerr << "Failed to write step " << t << std::endl;
      exit(2);
    }
  }
  return 0;
}
