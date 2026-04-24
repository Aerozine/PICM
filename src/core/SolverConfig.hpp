#pragma once
struct SolverConfig {
  enum class Type   { JACOBI, GAUSS_SEIDEL, RB_GS, MICCG0, CG };
  enum class Method { SL, VanillaPIC, FLIP, Mixed_FLIP_PIC, APIC };

  Method method    = Method::SL;
  Type   type      = Type::MICCG0;
  int    maxIters  = 1000;
  varType tolerance = 1e-4;
};
