#pragma once
#include "SceneObjects.hpp"
#include <nlohmann/json.hpp>
#include <ostream>
#include <string>

/**
 * @file Parameters.hpp
 * @brief Simulation configuration loaded from a JSON file.
 */

// Forward declaration — avoids pulling Fields2D into every translation unit
// that only needs grid dimensions or time-step values.
class Fields2D;

// SolverConfig
/**
 * @brief Configuration for the iterative pressure (Poisson) solver.
 */
struct SolverConfig {
  /// Available pressure solver algorithms.
  enum class Type {
    JACOBI,       ///< Jacobi iteration (parallelisable, slow convergence).
    GAUSS_SEIDEL, ///< Gauss-Seidel (faster convergence, sequential).
    RED_BLACK_GAUSS_SEIDEL ///< Red-black GS (parallelisable + fast
                           ///< convergence).
  };

  Type type = Type::GAUSS_SEIDEL; ///< Solver algorithm.
  int maxIters = 1000;            ///< Maximum number of iterations per step.
  double tolerance = 1e-2;        ///< Relative residual convergence threshold.

  /**
   * @brief Construct a SolverConfig from a JSON object.
   *
   * Recognised keys: @c "type", @c "max_iterations", @c "tolerance".
   * Unknown solver types fall back to GAUSS_SEIDEL with a warning.
   *
   * @param j JSON object node.
   * @return  Populated SolverConfig.
   */
  [[nodiscard]] static SolverConfig fromJson(const nlohmann::json &j);

  /// @return The solver type as a lowercase string (matches JSON key values).
  [[nodiscard]] std::string typeName() const;
};

// Parameters
/**
 * @brief All simulation parameters parsed from a JSON configuration file.
 *
 * ## Deferred scene construction
 * Scene objects (velocity patches, solid regions) are stored as raw JSON
 * subtrees and are **not** materialised into @c SceneObject instances until
 * @c applyToFields() is called. This keeps Parameters lightweight and avoids
 * any dependency on @c Fields2D in this header.
 */
class Parameters {
public:
  //Grid & time
  double dx = 0.01; ///< Cell width  in x (m).
  double dy = 0.01; ///< Cell height in y (m).
  double dt = 1e-4; ///< Time-step size (s).
  int nx = 100;     ///< Number of pressure cells in x.
  int ny = 100;     ///< Number of pressure cells in y.
  int nt = 100;     ///< Total number of time steps to simulate.

  // Physics
  double density = 1000.0; ///< Fluid density (kg/m³).

  // Output
  int sampling_rate = 1;          ///< Write output every N steps.
  std::string folder = "results"; ///< Output directory.
  std::string filename =
      "simulation"; ///< Base filename (unused at runtime, reserved).

  bool write_u = true;              ///< Write u-velocity field.
  bool write_v = true;              ///< Write v-velocity field.
  bool write_p = true;              ///< Write pressure field.
  bool write_div = false;           ///< Write divergence field (diagnostic).
  bool write_norm_velocity = false; ///< Write velocity magnitude (diagnostic).

  // Solver
  SolverConfig solver; ///< Pressure solver settings.

  // Life cycle
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
   * This is the only place where @c SceneObject instances are created.
   * Call once from the solver constructor after @c Fields2D is initialised.
   *
   * @param fields Target fields to mutate (velocities, solid labels).
   */
  void applyToFields(Fields2D &fields) const;

  /// Pretty-print all parameters to @p os (debug builds).
  friend std::ostream &operator<<(std::ostream &os, const Parameters &p);

private:
  // Raw JSON subtrees — SceneObjects are created lazily in applyToFields().
  nlohmann::json velocityU_json; ///< JSON node for initial u-velocity patches.
  nlohmann::json velocityV_json; ///< JSON node for initial v-velocity patches.
  nlohmann::json solid_json;     ///< JSON node for solid geometry.

  /**
   * @brief Populate members from a parsed JSON object.
   * @param j Root JSON object of the config file.
   */
  void loadFromJson(const nlohmann::json &j);

  /// Print command-line usage to stdout.
  static void printUsage(const char *prog);
};
