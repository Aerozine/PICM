#include "Parameters.hpp"
#include "Fields.hpp"
#include <cstring>
#include <fstream>
#include <iostream>

// SolverConfig
SolverConfig SolverConfig::fromJson(const nlohmann::json &j) {
  SolverConfig cfg;

  cfg.maxIters = static_cast<int>(j.value("max_iterations", 1000.0));
  cfg.tolerance = j.value("tolerance", 1e-4);
  cfg.type = Type::MICCG0;
  cfg.method = Method::SL;

  if (j.contains("type")) {
    const std::string t = j["type"].get<std::string>();
    if (t == "jacobi")
      cfg.type = Type::JACOBI;
    else if (t == "gauss_seidel")
      cfg.type = Type::GAUSS_SEIDEL;
    else if (t == "red_black_gauss_seidel")
      cfg.type = Type::RB_GS;
    else if (t == "miccg0")
      cfg.type = Type::MICCG0;
    else if (t == "cg")
      cfg.type = Type::CG;
    else
      std::cerr << "[SolverConfig] Unknown solver type '" << t
                << "' – defaulting to miccg0.\n";
  }

  if (j.contains("method"))
    cfg.method = methodFromJson(j["method"].get<std::string>());

  return cfg;
}

std::string SolverConfig::typeName() const {
  switch (type) {
  case Type::JACOBI:
    return "jacobi";
  case Type::GAUSS_SEIDEL:
    return "gauss_seidel";
  case Type::RB_GS:
    return "red_black_gauss_seidel";
  case Type::MICCG0:
    return "miccg0";
  case Type::CG:
    return "cg";
  }
  return "unknown";
}

SolverConfig::Method SolverConfig::methodFromJson(const nlohmann::json &j) {
  const std::string s = j.get<std::string>();
  if (s == "semilagrangian" || s == "sl")
    return Method::SL;
  if (s == "vanilla_pic" || s == "pic")
    return Method::VanillaPIC;
  if (s == "flip")
    return Method::FLIP;
  if (s == "mixed_flip_pic")
    return Method::Mixed_FLIP_PIC;
  if (s == "apic")
    return Method::APIC;
  std::cerr << "[SolverConfig] Unknown method '" << s
            << "' – defaulting to semi_lagrangian.\n";
  return Method::SL;
}

std::string SolverConfig::methodName(Method m) {
  switch (m) {
  case Method::SL:
    return "semi_lagrangian";
  case Method::VanillaPIC:
    return "vanilla_pic";
  case Method::FLIP:
    return "flip";
  case Method::Mixed_FLIP_PIC:
    return "mixed_flip_pic";
  case Method::APIC:
    return "apic";
  }
  return "unknown";
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
  // Output flags
  write_u = j.value("write_u", write_u);
  write_v = j.value("write_v", write_v);
  write_p = j.value("write_p", write_p);
  write_div = j.value("write_div", write_div);
  write_norm_velocity = j.value("write_norm_velocity", write_norm_velocity);
  write_smoke = j.value("write_smoke", write_smoke);
  write_countAliveParticles = j.value("write_countAliveParticles", write_countAliveParticles);
  ppcy = j.value("ppcy", ppcy);
  ppcx = j.value("ppcx", ppcx);

  freeSurface = j.value("freeSurface", freeSurface);
  gravity = j.value("gravity", gravity);
  refill = j.value("refill", refill);
  write_particles = j.value("write_particles", write_particles);
  // Output paths
  folder = j.value("folder", folder);
  // Boundary condition
  if (j.contains("velocityu"))
    velocityU_json = j["velocityu"];
  if (j.contains("velocityv"))
    velocityV_json = j["velocityv"];
  if (j.contains("pressure"))
    pressure_json = j["pressure"];
  if (j.contains("solid"))
    solid_json = j["solid"];
  if (j.contains("air"))
    air_json = j["air"];
  if (j.contains("fluid"))
    fluid_json = j["fluid"];
  if (j.contains("smoke"))
    smoke_json = j["smoke"];
  if (j.contains("solver")) {
    solver = SolverConfig::fromJson(j["solver"]);
    if (j.contains("method"))
      solver.method =
          SolverConfig::methodFromJson(j["method"].get<std::string>());
  }
}

void Parameters::applyToFields(Fields2D &fields) const {
  const std::map<std::string, int> vars = {{"nx", nx}, {"ny", ny}};

  if (!air_json.is_null()) {
    for (const auto &obj : parseSceneObjects(air_json, vars))
      obj->applyAir(fields);
  }
  if (!fluid_json.is_null()) {
    for (const auto &obj : parseSceneObjects(fluid_json, vars))
      obj->applyFluid(fields);
  }
  if (!smoke_json.is_null()) {
    for (const auto &obj : parseSceneObjects(smoke_json, vars))
      obj->applySmoke(fields);
  }
  if (!solid_json.is_null()) {
    for (const auto &obj : parseSceneObjects(solid_json, vars))
      obj->applySolid(fields);
  }

  if (!velocityU_json.is_null()) {
    for (const auto &obj : parseSceneObjects(velocityU_json, vars))
      obj->applyVelocityU(fields);
  }
  if (!velocityV_json.is_null()) {
    for (const auto &obj : parseSceneObjects(velocityV_json, vars))
      obj->applyVelocityV(fields);
  }
  if (!pressure_json.is_null()) {
    for (const auto &obj : parseSceneObjects(pressure_json, vars))
      obj->applyPressure(fields);
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
    DBG_PRINTF("[Parameters] Loaded from %s", path.c_str());
    return true;
  } catch (const std::exception &e) {
    std::cerr << "[Parameters] JSON parse error: " << e.what() << '\n';
    return false;
  }
}

bool Parameters::parseCommandLine(const int argc, char *argv[]) {
  if (argc == 2)
    return loadFromFile(argv[1]);
  printUsage(argv[0]);
  return false;
}

void Parameters::printUsage(const char *prog) {
  std::cout << "Usage: " << prog << " <config.json>\n";
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
