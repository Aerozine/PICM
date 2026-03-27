#include "core/Parameters.hpp"
#include "solvers/PIC/PIC.hpp"
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

#ifndef NDEBUG
    // Display parameters
    std::cout << params << std::endl;
#endif

    if (params.method == "pic") {
      PIC solver(params);
      solver.Run();
    } else {
      // Default: "semi_lagrangian"
      if (params.method != "semi_lagrangian")
        std::cerr << "[main] Unknown method '" << params.method
                  << "' – defaulting to semi_lagrangian.\n";
      SemiLagrangian solver(params);
      solver.Run();
    }

    std::cout << "Simulation completed successfully!\n";
    return 0;
  }
}