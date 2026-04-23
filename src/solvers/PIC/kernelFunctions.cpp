#include "PIC.hpp"
varType PIC::hat(varType r) const {
#ifdef USE_SPEED
    if (r >= varType(0) && r <= varType(1))
        return varType(1) - r;
    if (r >= varType(-1) && r < varType(0))
        return varType(1) + r;
    return varType(0);
#else
    if (varType(-1.5) <= r && r < varType(-0.5))
        return varType(0.5) * (r + varType(3.0 / 2.0)) * (r + varType(1.5));
    if (varType(-0.5) <= r && r < varType(0.5))
        return varType(0.75) - r * r;
    if (varType(0.5) <= r && r < varType(1.5))
        return varType(0.5) * (varType(1.5) - r) * (varType(1.5) - r);
    return varType(0);
#endif
}

varType PIC::dhat(varType r) const {
#ifdef USE_SPEED
    if (r > varType(0) && r < varType(1))
        return -varType(1);
    if (r > varType(-1) && r < varType(0))
      return varType(1);
    return varType(0);
#else
    if (varType(-1.5) <= r && r < varType(-0.5))
        return r + static_cast<varType>(1.5);
    if (varType(-0.5) <= r && r < varType(0.5))
        return -static_cast<varType>(2.0) * r;
    if (varType(0.5) <= r && r < varType(1.5))
        return r - static_cast<varType>(1.5);
    return static_cast<varType>(0);
#endif
}


