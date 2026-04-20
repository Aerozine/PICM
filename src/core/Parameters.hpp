#pragma once
#include "SceneObjects.hpp"
#include "SolverConfig_io.hpp"   // host-only: pulls in SolverConfig + json + string
#include <nlohmann/json.hpp>
#include <string>

/// @brief All simulation parameters parsed from a JSON configuration file.
class Parameters {
public:
  // Grid & time
  double dx = 0.01;
  double dy = 0.01;
  double dt = 1e-4;
  int    nx = 100;
  int    ny = 100;
  int    nt = 100;
  int ppcx  = 1;
  int ppcy  = 1;

  // Physics
  double density = 1000.0;

  // Output
  int         sampling_rate = 1;
  std::string folder        = "results";

  varType gravity = 0.0;

  bool write_u             = true;
  bool write_v             = true;
  bool write_p             = true;
  bool write_div           = false;
  bool write_norm_velocity = false;
  bool write_smoke         = false;
  bool write_particles     = false;
  bool refill              = false;

  /// true  → initialise the domain as AIR  (free-surface simulation)
  /// false → initialise as FLUID (fully-filled, incompressible default)
  bool freeSurface = false;

  SolverConfig solver;

  Parameters() = default;

  bool parseCommandLine(int argc, char *argv[]);
  bool loadFromFile(const std::string &path);
  void applyToFields(Fields2D &fields) const;

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

  void loadFromJson(const nlohmann::json &j);
  static void printUsage(const char *prog);
};
