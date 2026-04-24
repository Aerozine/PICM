#include "PIC.hpp"
// gravity through particles -> change order of functions in Step() in PIC.cpp
/*
void PIC::ApplyGravity() const {
  varType g = params.gravity;
OMP_PRAGMA(omp parallel for)
for (int idx = 0; idx < particles->size(); idx++) {
  varType v = particles->GetV(idx);
  particles->SetV(idx, v - dt * g);
}
}
// */