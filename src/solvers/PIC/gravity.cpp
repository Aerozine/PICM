#include "PIC.hpp"

void PIC::ApplyGravity() {
  for (int idx = 0; idx < particles->size(); ++idx) {
    varType g = 9.81;
    varType v = particles->GetV(idx);
    particles->SetV(idx, v - dt * g);
  }
}