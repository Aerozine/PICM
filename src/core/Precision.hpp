#pragma once
// Define in the cmake environment as option
#ifdef USE_FLOAT
typedef float Real;
#define REAL_EPSILON 1e-6f
// for compiler
#define REAL_LITERAL(x) x##f
#elif defined(USE_DOUBLE)
typedef double Real;
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
