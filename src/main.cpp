#include "core/Parameters.hpp"
#include "solvers/PIC/PIC.hpp"
#include "solvers/PIC/FLIP.hpp"
#include "solvers/PIC/APIC.hpp"
#include "solvers/SemiLagrangian/SemiLagrangian.hpp"

#include <iostream>

int main(int argc, char *argv[]) {
#ifndef NDEBUG
  std::cout << "Compiled with debug mode" << std::endl;
#endif

  // Parse parameters from command line
  Parameters params;
  if (!params.parseCommandLine(argc, argv)) {
    return 42;
  }
  std::cout << params << std::endl;

  // Initialize the solver
  if (params.method == "pic") {
    std::cout << "Particle In cell method \n" << std::endl;

    if (params.particleMethod == "vanilla_pic"){
      PIC solver(params);
      std::cout << "variant of the method: vanilla pic \n" << std::endl;
      solver.Run();
    } else if (params.particleMethod == "flip"){
      FLIP solver(params);
      std::cout << "variant of the method: flip \n" << std::endl;
      solver.Run();
    } else if (params.particleMethod == "apic"){
      APIC solver(params);
      std::cout << "variant of the method: apic \n" << std::endl;
      solver.Run();
    } else{
      std::cerr << "[main] Unknown PIC variant '" << params.particleMethod
                << "' – defaulting to vanilla pic.\n";
      PIC solver(params);
      solver.Run();
    }

  } else {
    // Default: "semi_lagrangian"
    if (params.method != "semi_lagrangian")
      std::cerr << "[main] Unknown method '" << params.method
                << "' – defaulting to semi_lagrangian.\n";

    std::cout << "SL" << std::endl;
    SemiLagrangian solver(params);
    solver.Run();
  }
  std::cout << "Simulation completed successfully!\n";
  return 0;
}