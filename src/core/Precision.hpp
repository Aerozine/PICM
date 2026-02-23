#pragma once
#include <omp.h>

/**
 * @file Precision.hpp
 * @brief Compile-time floating-point precision selection.
 *
 * Define either USE_FLOAT or USE_DOUBLE via CMake to select the working
 * precision for the entire simulation. All numerical fields, grids, and
 * solver variables use @c varType.
 */

#ifdef USE_FLOAT

using varType = float; ///< Simulation floating-point type (32-bit).

#define REAL_EPSILON 1e-6f   ///< Small epsilon for float comparisons.
#define REAL_LITERAL(x) x##f ///< Suffix literal with 'f' for float precision.
#define PRECISION_STRING "float (32-bit)"

#elif defined(USE_DOUBLE)

using varType = double; ///< Simulation floating-point type (64-bit).

#define REAL_EPSILON 1e-15 ///< Small epsilon for double comparisons.
#define REAL_LITERAL(x) x  ///< No suffix needed for double precision.
#define PRECISION_STRING "double (64-bit)"

#else
#error "Precision not defined: set either USE_FLOAT or USE_DOUBLE in CMake."
#endif

/// @brief Wall-clock time in seconds (via OpenMP).
#define GET_TIME() (omp_get_wtime())
