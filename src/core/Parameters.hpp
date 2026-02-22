#pragma once
#include "SceneObjects.hpp"
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

struct SolverConfig {
  // stores more efficiently solver type instead of strcmp string name(may be
  // critical due to the check in the loop)
  enum class Type { JACOBI, GAUSS_SEIDEL };
  // Default values
  Type type = Type::GAUSS_SEIDEL;
  int maxIters = 1000;
  double tolerance = 1e-2;

  static SolverConfig fromJson(const nlohmann::json &j);
  std::string typeName() const;
};
// TODO maybe not store Velocity and Solid config INSIDE
// the params structure
// maybe conserving the file inside params and when loading read the file ?
struct VelocityConfig {
  std::vector<std::unique_ptr<SceneObject>> objects;

  static VelocityConfig fromJson(const nlohmann::json &j,
                                 const std::map<std::string, int> &vars);
  VelocityConfig() = default;
  VelocityConfig(VelocityConfig &&) = default;
  VelocityConfig &operator=(VelocityConfig &&) = default;
};

//  Solid block
struct SolidConfig {
  std::vector<std::unique_ptr<SceneObject>> objects;

  static SolidConfig fromJson(const nlohmann::json &j,
                              const std::map<std::string, int> &vars);
  SolidConfig() = default;
  SolidConfig(SolidConfig &&) = default;
  SolidConfig &operator=(SolidConfig &&) = default;
};

//  Initial condition type
//
//  "custom"       – use velocityu / velocityv / solid JSON blocks (default)
//  "taylor_green" – WIP  call Fields2D::InitTaylorGreen(); ignore velocity
//  blocks
//                   and solid blocks (periodic domain, no walls)
enum class InitialCondition { CUSTOM, TAYLOR_GREEN };

InitialCondition parseInitialCondition(const std::string &s);

//  Top-level parameters
class Parameters {
public:
  // default values
  double dx = 0.01, dy = 0.01, dt = 1e-12;
  int nx = 100, ny = 100, nt = 100;
  int sampling_rate = 1;
  double density = 1000;

  // flags
  bool write_u = true;
  bool write_v = true;
  bool write_p = true;
  bool write_div = false;
  bool write_norm_velocity = false;
  
  // Source condition
  bool entryFlow = false;
  double entryFlowVelocity = 0.0;

  // ── Output paths ─────────────────────────────────────────────────────────
  std::string folder = "results";
  std::string filename = "simulation";

  // ── Initial condition ─────────────────────────────────────────────────────
  // Selects which init path the solver takes.
  InitialCondition initialCondition = InitialCondition::CUSTOM;

  // Optional amplitude for taylor_green (default 1.0)
  double taylorGreenAmplitude = 1.0;

  // ── Scene / solver config (used when initialCondition == CUSTOM) ──────────
  VelocityConfig velocityU;
  VelocityConfig velocityV;
  SolidConfig solid;
  SolverConfig solver;

  // ── I/O ──────────────────────────────────────────────────────────────────
  Parameters();
  bool loadFromFile(const std::string &path);
  bool parseCommandLine(int argc, char *argv[]);

  friend std::ostream &operator<<(std::ostream &os, const Parameters &p);

private:
  void setDefaults();
  void loadFromJson(const nlohmann::json &j);
  static void printUsage(const char *prog);
};
