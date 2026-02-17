#include "core/Parameters.hpp"
#include "solvers/SemiLagrangian/SemiLagrangian.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
#ifndef NDEBUG
  std::cout << "Compiled with debug mode" << std::endl;
#endif

  // Parse parameters from command line
  Parameters params;
  if (!params.parseCommandLine(argc, argv)) {
    return 1;
  }

#ifndef NDEBUG
  // Display parameters
  std::cout << params << std::endl;
#endif

  // Create and run solver
  SemiLagrangian solver(params);
  solver.Run();

  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}
