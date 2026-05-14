#pragma once

#include "Precision.hpp"
#include <nlohmann/json.hpp>
#include <string>

struct SolverConfig {
  enum class Type { JACOBI, GAUSS_SEIDEL, RB_GS, MICCG0, CG };
  enum class Method {
    SL,
    VanillaPIC,
    FLIP,
    Mixed_FLIP_PIC,
    APIC
  };

  Method method = Method::SL;
  Type type = Type::MICCG0;
  bool surfaceTension = false;
  int maxIters = 1000;
  varType tolerance = 1e-4;
};

SolverConfig solverConfigFromJson(const nlohmann::json &j);
std::string solverTypeName(const SolverConfig &cfg);
SolverConfig::Method solverMethodFromJson(const nlohmann::json &j);
std::string solverMethodName(SolverConfig::Method m);
