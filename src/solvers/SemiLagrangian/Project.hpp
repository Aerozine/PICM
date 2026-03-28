#pragma once
#include "../../core/Parameters.hpp"
#include "../../core/IterativeMethods.hpp"
void MakeIncompressible(const Parameters & p,Fields2D & f );

void updateVelocities(const Parameters & p,Fields2D & f );

void solvePressure(const Parameters & p,Fields2D & f );
