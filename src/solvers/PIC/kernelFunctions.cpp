#include "PIC.hpp"

varType PIC::hat(varType r) const {
    if (params.kernelOrder == 1)
        return hatTrilinear(r);
    else if (params.kernelOrder == 2)
        return hatQuadraticBSpline(r);
    else
        return varType(0);
}

varType PIC::dhat(varType r) const {
    if (params.kernelOrder == 1)
        return dhatTrilinear(r);
    else if (params.kernelOrder == 2)
        return dhatQuadraticBSpline(r);
    else
        return varType(0);
}

// @todo handle different hat correctly
// needs to be handle in the Cmake for precision over cost
varType PIC::hatTrilinear(varType r) const {
  if (r >= varType(0) && r <= varType(1))
    return varType(1) - r;
  else if (r >= varType(-1) && r < varType(0))
    return varType(1) + r;
  else
    return varType(0);
}

varType PIC::dhatTrilinear(varType r) const {
    if (r > varType(0) && r < varType(1))
        return -varType(1);
    else if (r > varType(-1) && r < varType(0))
        return varType(1);
    else
        return varType(0);
}

varType PIC::hatQuadraticBSpline(varType r) const {
  if (varType(-1.5) <= r && r < varType(-0.5))
    return varType(0.5) * (r + varType(3.0 / 2.0)) * (r + varType(3.0 / 2.0));
  if (varType(-0.5) <= r && r < varType(0.5))
    return varType(0.75) - r * r;
  if (varType(0.5) <= r && r < varType(1.5))
    return varType(0.5) * (varType(3.0 / 2.0) - r) * (varType(3.0 / 2.0) - r);
  return varType(0);
}

varType PIC::dhatQuadraticBSpline(varType r) const {
  if (varType(-1.5) <= r && r < varType(-0.5))
    return r + static_cast<varType>(1.5);
  if (varType(-0.5) <= r && r < varType(0.5))
    return -static_cast<varType>(2.0) * r;
  if (varType(0.5) <= r && r < varType(1.5))
    return r - static_cast<varType>(1.5);
  return static_cast<varType>(0);
}

