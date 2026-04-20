#pragma once

#include "SolverConfig.hpp"
#include <nlohmann/json.hpp>
#include <string>

SolverConfig   solverConfigFromJson(const nlohmann::json &j);
std::string    solverTypeName(const SolverConfig &cfg);
SolverConfig::Method solverMethodFromJson(const nlohmann::json &j);
std::string    solverMethodName(SolverConfig::Method m);
