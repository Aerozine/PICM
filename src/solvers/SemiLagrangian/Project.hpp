#pragma once
#include "../../core/IterativeMethods.hpp"
#include "../../core/Parameters.hpp"
void MakeIncompressible(const Parameters &p, Fields2D &f);

void updateVelocities(const Parameters &p, Fields2D &f);

void solvePressure(const Parameters &p, Fields2D &f);
