#include "core/Parameters.hpp"
#include "solvers/APIC/APIC.hpp"
#include "solvers/FLIP/FLIP.hpp"
#include "solvers/PIC/PIC.hpp"
#include "solvers/SemiLagrangian/SemiLagrangian.hpp"
#ifdef USE_CUDA
#include "solvers/GPIC/GPIC.hpp"
#endif

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
  switch (params.solver.method) {
  case SolverConfig::Method::SL: {
    std::cout << "SL\n";
    SemiLagrangian solver(params);
    solver.Run();
    break;
  }
  case SolverConfig::Method::VanillaPIC: {
    std::cout << "VanillaPIC\n";
    PIC solver(params);
    solver.Run();
    break;
  }
  case SolverConfig::Method::FLIP: {
    std::cout << "FLIP\n";
    FLIP solver(params);
    solver.Run();
    break;
  }
  case SolverConfig::Method::Mixed_FLIP_PIC:
    std::cerr << "Mixed_FLIP_PIC not yet implemented\n";
    return 1;
  case SolverConfig::Method::APIC: {
    std::cout << "APIC\n";
    APIC solver(params);
    solver.Run();
    break;
  }
#ifdef USE_CUDA
  case SolverConfig::Method::GPIC:
  case SolverConfig::Method::GFLIP:
  case SolverConfig::Method::GAPIC: {
    std::cout << solverMethodName(params.solver.method) << " (GPU)\n";
    GPIC solver(params);
    solver.Run();
    break;
  }
#else
  case SolverConfig::Method::GPIC:
  case SolverConfig::Method::GFLIP:
  case SolverConfig::Method::GAPIC:
    std::cerr << "[" << solverMethodName(params.solver.method)
              << "] Requires -DUSE_GPU=ON build.\n";
    return 1;
#endif
  }
  std::cout << "Simulation completed successfully!\n";
  return 0;
}
