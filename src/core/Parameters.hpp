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
    MICCG0        ///< Modified Incomplete Cholesky CG level 0 ( better )
  };

  Type type;
  int maxIters;     ///< Maximum number of iterations per step.
  double tolerance; ///< Relative residual convergence threshold.

  /**
   * @brief Construct a SolverConfig from a JSON object.
   *
   * @param j JSON object node.
   * @return  Populated SolverConfig.
   */
  static SolverConfig fromJson(const nlohmann::json &j);

  /// @return The solver type as a lowercase string (matches JSON key values).
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
  int ppcx = 1;
  int ppcy = 1;

  // Physics
  double density = 1000.0; ///< Fluid density (kg/m3).

  // Output
  int sampling_rate = 1;          ///< Write output every N steps.
  std::string folder = "results"; ///< Output directory.
  std::string filename =
      "simulation"; ///< Base filename (unused at runtime, reserved).

  varType gravity = 0.0;

  bool write_u = true;
  bool write_v = true;
  bool write_p = true;
  bool write_div = false;
  bool write_norm_velocity = false;
  bool write_smoke = false;
  bool write_particles = false;
  bool refill = false;

  std::string method = "semi_lagrangian";
  std::string freeSurface = "no";
  std::string particleMethod = "vanilla_pic";
  SolverConfig solver;

  Parameters() = default;

  /**
   * @brief Parse @c -c / @c --config \<path\> from @c argv and load the file.
   * @param argc Argument count from @c main.
   * @param argv Argument vector from @c main.
   * @return @c true on success, @c false on error (usage is printed).
   */
  bool parseCommandLine(int argc, char *argv[]);

  /**
   * @brief Load parameters from a JSON file.
   * @param path Path to the .json config file.
   * @return @c true on success, @c false if the file could not be opened or
   *         parsed.
   */
  bool loadFromFile(const std::string &path);

  /**
   * @brief Instantiate scene objects from the stored JSON and apply them to
   *        @p fields, then immediately discard the temporary objects.
   *
   * @param fields Target fields to mutate (velocities, solid labels).
   */
  void applyToFields(Fields2D &fields) const;

  /// Pretty-print all parameters to @p os (debug builds).
  friend std::ostream &operator<<(std::ostream &os, const Parameters &p);

private:
  nlohmann::json velocityU_json;
  nlohmann::json velocityV_json;
  nlohmann::json pressure_json;
  nlohmann::json solid_json;
  nlohmann::json air_json;
  nlohmann::json fluid_json;
  nlohmann::json smoke_json;
  nlohmann::json particles_json;

  /**
   * @brief Populate members from a parsed JSON object.
   * @param j Root JSON object of the config file.
   */
  void loadFromJson(const nlohmann::json &j);

  /// Print command-line usage to stdout.
  static void printUsage(const char *prog);
};
