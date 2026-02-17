#pragma once
#include <omp.h>
// Define in the cmake environment as option
#ifdef USE_FLOAT
// TODO remove real of no other use
typedef float Real;
typedef float varType;

// ?
#define REAL_EPSILON 1e-6f
// for compiler
// TODO remove also
#define REAL_LITERAL(x) x##f
#elif defined(USE_DOUBLE)
typedef double Real;
typedef double varType;
// for compiler
#define REAL_LITERAL(x) x
#else
#error "Must define either USE_FLOAT or USE_DOUBLE"
#endif

// Optional: Print which precision is being used
#ifdef USE_FLOAT
#define PRECISION_STRING "float (32-bit)"
#else
#define PRECISION_STRING "double (64-bit)"
#endif

// #define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time
#define GET_TIME() (omp_get_wtime()) // wall time