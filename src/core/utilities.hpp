#pragma once
#include "Grid2D.hpp"

typedef float varType

Grid2D Div(Grid2D field);

varType InterpolateOneValue(Grid2D field, size_t x, size_t y);

