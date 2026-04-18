#pragma once
#include <nlohmann/json.hpp>
#include <string>

struct SolverConfig {
  enum class Type { JACOBI, GAUSS_SEIDEL, RB_GS, MICCG0, CG };
  enum class Method { SL, VanillaPIC, FLIP, Mixed_FLIP_PIC, APIC };

  static Method methodFromJson(const nlohmann::json &j);
  static std::string methodName(Method m);
  static SolverConfig fromJson(const nlohmann::json &j);
  std::string typeName() const;

  Method method = Method::SL;
  Type type = Type::MICCG0;
  int maxIters = 1000;
  double tolerance = 1e-4;
};
