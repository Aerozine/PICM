#pragma once
#include "../../core/Parameters.hpp"
#include "../Solver.hpp"

class SemiLagrangian : public Solver {
public:
  /**
   * @brief Construct the solver, initialise fields, and open output writers.
   * @param params Simulation parameters (non-owning reference, must outlive
   *               this object).
   */

  explicit SemiLagrangian(Parameters &params);

  ~SemiLagrangian() override = default;

  // SemiLagrangian(const SemiLagrangian &) = delete;
  // SemiLagrangian &operator=(const SemiLagrangian &) = delete;

  /// @brief Run the full simulation loop (nt steps) and write output.
  void Run() override;

  /// @brief Advance the simulation by one time step.
  void Step() override;
};
