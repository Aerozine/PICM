#include "Parameters.hpp"
#include "Fields.hpp"
#include <cstring>
#include <fstream>
#include <iostream>
#include <string_view>

// SolverConfig
SolverConfig SolverConfig::fromJson(const nlohmann::json &j) {
    SolverConfig cfg{
        .type      = Type::RED_BLACK_GAUSS_SEIDEL,
        .maxIters  = j.value("max_iterations", 1000),
        .tolerance = j.value("tolerance", 1e-4),
    };

    if (j.contains("type")) {

        const std::string t = j["type"].get<std::string>();

        if (t == "jacobi") cfg.type = Type::JACOBI;
        else if (t == "gauss_seidel") cfg.type = Type::GAUSS_SEIDEL;
        else if (t == "red_black_gauss_seidel") cfg.type = Type::RED_BLACK_GAUSS_SEIDEL;
        else
            std::cerr << "[SolverConfig] Unknown solver type '" << t
                      << "' – defaulting to red_black_gauss_seidel.\n";
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

void Parameters::loadFromJson(const nlohmann::json &j) {
  // dx is set with the value dx found in the json or keep its previous value
  dx = j.value("dx", dx);
  dy = j.value("dy", dy);
  dt = j.value("dt", dt);
  nx = j.value("nx", nx);
  ny = j.value("ny", ny);
  nt = j.value("nt", nt);
  sampling_rate = j.value("sampling_rate", sampling_rate);
  density = j.value("density", density);
  source = j.value("source", source);
  // Output flags
  write_u = j.value("write_u", write_u);
  write_v = j.value("write_v", write_v);
  write_p = j.value("write_p", write_p);
  write_div = j.value("write_div", write_div);
  write_norm_velocity = j.value("write_norm_velocity", write_norm_velocity);
  write_smoke = j.value("write_smoke", write_smoke);
  // Output paths
  folder = j.value("folder", folder);
  filename = j.value("filename", filename);
  // Boundary condition
  if (j.contains("velocityu"))
    velocityU_json = j["velocityu"];
  if (j.contains("velocityv"))
    velocityV_json = j["velocityv"];
  if (j.contains("solid"))
    solid_json = j["solid"];
  if (j.contains("smoke"))
    smoke_json = j["smoke"];
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
DBG_PRINTF("[Parameters] Loaded from %s",path.c_str());
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
// allow us to use the << operator to print params
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
