#include "Parameters.hpp"
#include "Fields.hpp"
#include <cstring>
#include <fstream>
#include <iostream>
#include <string_view>

// SolverConfig

SolverConfig SolverConfig::fromJson(const nlohmann::json &j) {
  SolverConfig cfg;

  if (j.contains("max_iterations"))
    cfg.maxIters = j["max_iterations"].get<int>();

  if (j.contains("tolerance"))
    cfg.tolerance = j["tolerance"].get<double>();

  if (j.contains("type")) {
    const std::string t = j["type"].get<std::string>();
    if (t == "jacobi")
      cfg.type = Type::JACOBI;
    else if (t == "gauss_seidel")
      cfg.type = Type::GAUSS_SEIDEL;
    else if (t == "red_black_gauss_seidel")
      cfg.type = Type::RED_BLACK_GAUSS_SEIDEL;
    else
      std::cerr << "[SolverConfig] Unknown solver type '" << t
                << "' – defaulting to gauss_seidel.\n";
  }
  return cfg;
}

std::string SolverConfig::typeName() const {
  switch (type) {
  case Type::JACOBI:
    return "jacobi";
  case Type::GAUSS_SEIDEL:
    return "gauss_seidel";
  case Type::RED_BLACK_GAUSS_SEIDEL:
    return "red_black_gauss_seidel";
  }
  return "unknown"; // unreachable, silences -Wreturn-type
}

// Parameters

void Parameters::loadFromJson(const nlohmann::json &j) {
  // Helper lambda: assign a field only if the key is present in the JSON.
  // Using a lambda avoids repeating the j.contains / j[key].get<T>() pattern.
  auto load = [&j](const char *key, auto &member) {
    if (j.contains(key))
      member = j[key].get<std::decay_t<decltype(member)>>();
  };

  // Grid & time
  load("dx", dx);
  load("dy", dy);
  load("dt", dt);
  load("nx", nx);
  load("ny", ny);
  load("nt", nt);
  load("sampling_rate", sampling_rate);
  load("density", density);

  // simulation condition
  load("source", source);

  // Output flags
  load("write_u", write_u);
  load("write_v", write_v);
  load("write_p", write_p);
  load("write_div", write_div);
  load("write_norm_velocity", write_norm_velocity);
  load("write_smoke", write_smoke);

  // Output paths
  load("folder", folder);
  load("filename", filename);

  // Scene geometry — store raw JSON; SceneObjects are built lazily in
  // applyToFields() so that Parameters has no dependency on Fields2D.
  if (j.contains("velocityu"))
    velocityU_json = j["velocityu"];
  if (j.contains("velocityv"))
    velocityV_json = j["velocityv"];
  if (j.contains("solid"))
    solid_json = j["solid"];
  if (j.contains("smoke"))
    smoke_json = j["smoke"];

  // Solver
  if (j.contains("solver"))
    solver = SolverConfig::fromJson(j["solver"]);
}

void Parameters::applyToFields(Fields2D &fields) const {
  const std::map<std::string, int> vars = {{"nx", nx}, {"ny", ny}};

  if (!velocityU_json.is_null()) {
    for (const auto &obj : parseSceneObjects(velocityU_json, vars))
      obj->applyVelocityU(fields);
  }
  if (!velocityV_json.is_null()) {
    for (const auto &obj : parseSceneObjects(velocityV_json, vars))
      obj->applyVelocityV(fields);
  }
  if (!solid_json.is_null()) {
    for (const auto &obj : parseSceneObjects(solid_json, vars))
      obj->applySolid(fields);
  }
  if (!smoke_json.is_null()) {
    for (const auto &obj : parseSceneObjects(smoke_json, vars))
      obj->applySmoke(fields);
  }
}

bool Parameters::loadFromFile(const std::string &path) {
  try {
    std::ifstream file(path);
    if (!file.is_open()) {
      std::cerr << "[Parameters] Could not open '" << path << "'\n";
      return false;
    }
    nlohmann::json j;
    file >> j;
    loadFromJson(j);
#ifndef NDEBUG
    std::cout << "[Parameters] Loaded from '" << path << "'\n";
#endif
    return true;
  } catch (const std::exception &e) {
    std::cerr << "[Parameters] JSON parse error: " << e.what() << '\n';
    return false;
  }
}

bool Parameters::parseCommandLine(int argc, char *argv[]) {
  // Expect exactly:  <prog> -c <path>  or  <prog> --config <path>
  if (argc == 3) {
    const std::string_view flag = argv[1];
    if (flag == "-c" || flag == "--config")
      return loadFromFile(argv[2]);
  }
  printUsage(argv[0]);
  return false;
}

void Parameters::printUsage(const char *prog) {
  // RTFM
  std::cout << "Usage: " << prog << " -c <config.json>\n";
}

std::ostream &operator<<(std::ostream &os, const Parameters &p) {
  os << "\n=== Simulation Parameters ===\n"
     << "  Grid    : " << p.nx << " x " << p.ny << "  dx=" << p.dx
     << "  dy=" << p.dy << '\n'
     << "  Time    : nt=" << p.nt << "  dt=" << p.dt << '\n'
     << "  Density : " << p.density << '\n'
     << "  Sampling: every " << p.sampling_rate << " step(s)" << '\n'
     << "  Solver  : " << p.solver.typeName()
     << "  maxIter=" << p.solver.maxIters << "  tol=" << p.solver.tolerance
     << '\n'
     << "  Output  : folder='" << p.folder << "'\n"
     << "  Write   : u=" << p.write_u << " v=" << p.write_v
     << " p=" << p.write_p << " div=" << p.write_div
     << " norm=" << p.write_norm_velocity << '\n'
     << "  InitVelU: " << (!p.velocityU_json.is_null() ? "defined" : "none")
     << '\n'
     << "  InitVelV: " << (!p.velocityV_json.is_null() ? "defined" : "none")
     << '\n'
     << "  smoke: " << (!p.smoke_json.is_null() ? "defined" : "none") << '\n'
     << "  Solid   : " << (!p.solid_json.is_null() ? "defined" : "none") << '\n'
     << "=============================\n";
  return os;
}
