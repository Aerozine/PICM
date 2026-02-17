#include "Parameters.hpp"
#include <cstring>
#include <fstream>
#include <iostream>

// specify differents built in IC by default ret custom
InitialCondition parseInitialCondition(const std::string &s) {
  if (s == "taylor_green")
    return InitialCondition::TAYLOR_GREEN;
  if (s == "custom")
    return InitialCondition::CUSTOM;
  std::cerr << "[Parameters] Unknown initialcondition '" << s
            << "' – defaulting to 'custom'.\n";
  return InitialCondition::CUSTOM;
}

SolverConfig SolverConfig::fromJson(const nlohmann::json &j) {
  SolverConfig cfg;
  if (j.contains("max_iterations"))
    cfg.maxIters = j["max_iterations"].get<int>();
  if (j.contains("tolerance"))
    cfg.tolerance = j["tolerance"].get<double>();
  if (j.contains("type")) {
    std::string t = j["type"].get<std::string>();
    if (t == "jacobi")
      cfg.type = Type::JACOBI;
    else if (t == "gauss_seidel")
      cfg.type = Type::GAUSS_SEIDEL;
    else
      std::cerr << "[Parameters] Unknown solver type '" << t
                << "' – defaulting to gauss_seidel.\n";
  }
  return cfg;
}
// just for printing purpose
std::string SolverConfig::typeName() const {
  switch (type) {
  case Type::JACOBI:
    return "jacobi";
  case Type::GAUSS_SEIDEL:
    return "gauss_seidel";
  default:
    return "unknown";
  }
}

VelocityConfig
VelocityConfig::fromJson(const nlohmann::json &j,
                         const std::map<std::string, int> &vars) {
  VelocityConfig cfg;
  cfg.objects = parseSceneObjects(j, vars);
  return cfg;
}

SolidConfig SolidConfig::fromJson(const nlohmann::json &j,
                                  const std::map<std::string, int> &vars) {
  SolidConfig cfg;
  cfg.objects = parseSceneObjects(j, vars);
  return cfg;
}

Parameters::Parameters() {}

void Parameters::loadFromJson(const nlohmann::json &j) {
  // scalar values
  if (j.contains("dx"))
    dx = j["dx"];
  if (j.contains("dy"))
    dy = j["dy"];
  if (j.contains("dt"))
    dt = j["dt"];
  if (j.contains("nx"))
    nx = j["nx"];
  if (j.contains("ny"))
    ny = j["ny"];
  if (j.contains("nt"))
    nt = j["nt"];
  if (j.contains("sampling_rate"))
    sampling_rate = j["sampling_rate"];
  if (j.contains("density"))
    density = j["density"];
  // writing flags
  if (j.contains("write_u"))
    write_u = j["write_u"];
  if (j.contains("write_v"))
    write_v = j["write_v"];
  if (j.contains("write_p"))
    write_p = j["write_p"];
  if (j.contains("write_div"))
    write_div = j["write_div"];
  if (j.contains("write_norm_velocity"))
    write_norm_velocity = j["write_norm_velocity"];
  // writing folder
  if (j.contains("folder"))
    folder = j["folder"].get<std::string>();
  if (j.contains("filename"))
    filename = j["filename"].get<std::string>();

  // Init condition
  if (j.contains("initialcondition")) {
    const auto &ic = j["initialcondition"];
    if (ic.contains("type"))
      initialCondition = parseInitialCondition(ic["type"].get<std::string>());
    if (ic.contains("amplitude"))
      taylorGreenAmplitude = ic["amplitude"].get<double>();
  }

  // build map of symbols
  // This allow us to say for example a rectangle located at nx/2
  const std::map<std::string, int> vars = {{"nx", nx}, {"ny", ny}};
  // load initial data
  if (j.contains("velocityu"))
    velocityU = VelocityConfig::fromJson(j["velocityu"], vars);
  if (j.contains("velocityv"))
    velocityV = VelocityConfig::fromJson(j["velocityv"], vars);
  if (j.contains("solid"))
    solid = SolidConfig::fromJson(j["solid"], vars);
  if (j.contains("solver"))
    solver = SolverConfig::fromJson(j["solver"]);
}

bool Parameters::loadFromFile(const std::string &path) {
  try {
    std::ifstream file(path);
    if (!file.is_open()) {
      std::cerr << "Error: Could not open '" << path << "'\n";
      return false;
    }
    nlohmann::json j;
    file >> j;
    loadFromJson(j);
#ifndef NDEBUG
    std::cout << "Loaded parameters from '" << path << "'\n";
#endif
    return true;
  } catch (const std::exception &e) {
    std::cerr << "Error parsing JSON: " << e.what() << '\n';
    return false;
  }
}

bool Parameters::parseCommandLine(int argc, char *argv[]) {
  if (argc == 1) {
    printUsage(argv[0]);
    return false;
  }
  for (int i = 1; i < argc; ++i) {
    if ((std::strcmp(argv[i], "-c") == 0 ||
         std::strcmp(argv[i], "--config") == 0) &&
        i + 1 < argc) {
      return loadFromFile(argv[++i]);
    } else {
      printUsage(argv[0]);
      return false;
    }
  }
  return true;
}

void Parameters::printUsage(const char *prog) {
  std::cout << "Usage: " << prog << " -c <config.json>\n";
}

std::ostream &operator<<(std::ostream &os, const Parameters &p) {
  std::string icName = (p.initialCondition == InitialCondition::TAYLOR_GREEN)
                           ? "taylor_green"
                           : "custom";
  os << "\n=== Simulation Parameters ===\n"
     << "  Grid    : " << p.nx << " x " << p.ny << "  dx=" << p.dx
     << " dy=" << p.dy << '\n'
     << "  Time    : nt=" << p.nt << " dt=" << p.dt << '\n'
     << "  Density : " << p.density << '\n'
     << "  Sampling: " << p.sampling_rate << '\n'
     << "  IC      : " << icName;
  if (p.initialCondition == InitialCondition::TAYLOR_GREEN)
    os << "  amplitude=" << p.taylorGreenAmplitude;
  os << '\n'
     << "  Solver  : " << p.solver.typeName()
     << "  maxIter=" << p.solver.maxIters << "  tol=" << p.solver.tolerance
     << '\n'
     << "  Output  : folder='" << p.folder << "'\n"
     << "  Write   : u=" << p.write_u << " v=" << p.write_v
     << " p=" << p.write_p << " div=" << p.write_div
     << " norm=" << p.write_norm_velocity << '\n'
     << "  velocityU objects : " << p.velocityU.objects.size() << '\n'
     << "  velocityV objects : " << p.velocityV.objects.size() << '\n'
     << "  solid    objects  : " << p.solid.objects.size() << '\n'
     << "=============================\n";
  return os;
}
