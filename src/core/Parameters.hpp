#pragma once
#include "SceneObjects.hpp"
#include <nlohmann/json.hpp>
#include <ostream>
#include <string>

class Fields2D;

/**
 * @brief Configuration for the iterative pressure (Poisson) solver.
 */
struct SolverConfig {
  /// Available pressure solver algorithms.
  enum class Type {
    JACOBI,       ///< Jacobi iteration (parallelisable, slow convergence).
    GAUSS_SEIDEL, ///< Gauss-Seidel (faster convergence, sequential).
    RB_GS,        ///< Red-black GS (parallelisable + fast convergence).
    MICCG0        ///< Modified Incomplete Cholesky CG level 0 (best).
  };

  Type type       = Type::MICCG0;
  int maxIters    = 1000;
  double tolerance = 1e-4;

  static SolverConfig fromJson(const nlohmann::json &j);

  /// @return The solver type as a lowercase string matching JSON key values.
  std::string typeName() const;
};

/// @brief All simulation parameters parsed from a JSON configuration file.
class Parameters {
public:
  // Grid & time
  double dx = 0.01; ///< Cell width  in x (m).
  double dy = 0.01; ///< Cell height in y (m).
  double dt = 1e-4; ///< Time-step size (s).
  int nx = 100;     ///< Number of pressure cells in x.
  int ny = 100;     ///< Number of pressure cells in y.
  int nt = 100;     ///< Total number of time steps to simulate.

  // PIC particle seeding (ignored by SemiLagrangian)
  int ppcx = 1; ///< Particles per cell in x.
  int ppcy = 1; ///< Particles per cell in y.

  // Physics
  double density = 1000.0; ///< Fluid density (kg/m3).

  // Output
  int sampling_rate = 1;            ///< Write output every N steps.
  std::string folder   = "results"; ///< Output directory.
  std::string filename = "simulation"; ///< Base filename (reserved).

  bool source = false; ///< Whether to use a source term.

  bool write_u             = true;
  bool write_v             = true;
  bool write_p             = true;
  bool write_div           = false;
  bool write_norm_velocity = false;
  bool write_smoke         = false;
  bool write_particles     = false;

  /// Solver selection: "semi_lagrangian" (default) or "pic".
  std::string method = "semi_lagrangian";

  SolverConfig solver;

  Parameters() = default;

  /**
   * @brief Parse @c -c / @c --config \<path\> from @c argv and load the file.
   */
  bool parseCommandLine(int argc, char *argv[]);

  /**
   * @brief Load parameters from a JSON file.
   */
  bool loadFromFile(const std::string &path);

  /**
   * @brief Instantiate scene objects from the stored JSON and apply them to
   *        @p fields.
   */
  void applyToFields(Fields2D &fields) const;

  friend std::ostream &operator<<(std::ostream &os, const Parameters &p);

private:
  nlohmann::json velocityU_json;
  nlohmann::json velocityV_json;
  nlohmann::json pressure_json;
  nlohmann::json solid_json;
  nlohmann::json smoke_json;

  void loadFromJson(const nlohmann::json &j);
  static void printUsage(const char *prog);
};
