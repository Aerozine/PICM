#pragma once
#include "SceneObjects.hpp"
#include "SolverConfig.hpp"
#include <nlohmann/json.hpp>
#include <string>

/// @brief All simulation parameters parsed from a JSON configuration file.
class Parameters {
public:
  // Grid & time
  varType dx = 0.01;
  varType dy = 0.01;
  varType dt = 1e-4;
  int nx = 100;
  int ny = 100;
  int nt = 100;
  int ppcx = 1;
  int ppcy = 1;

  // Physics
  varType density = 1000.0;

  // Output
  int sampling_rate = 1; ///< Write output every N steps.
  int kernelOrder = 2;
  std::string folder = "results"; ///< Output directory.

  varType gravity = 0.0;
  varType coefPic = 0.05;
  varType particleRadius = dx/10.0;
  varType gamma = 0.073;
  varType max_cfl = REAL_LITERAL(0.95);
  varType interactionStiffness = 0.0;
  varType interactionExponent = 0.0;

  bool write_u = true;
  bool write_v = true;
  bool write_p = true;
  bool write_div = false;
  bool write_norm_velocity = false;
  bool write_smoke = false;
  bool write_particles = false;
  bool refill = false;

  /// true  → initialise the domain as AIR  (free-surface simulation)
  /// false → initialise as FLUID (fully-filled, incompressible default)
  bool freeSurface = false;
  bool surfaceTension = false;

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
  void applySmoke(Grid2D &smoke, Fields2D &f);

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
