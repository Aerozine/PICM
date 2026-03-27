#pragma once
#include "Precision.hpp"
#include <vector>

/**
 * @file Particles.hpp
 * @brief Flat particle pool for the PIC solver.
 */

struct Vec2 {
  varType x = 0.0;
  varType y = 0.0;
};

struct Particle {
  Vec2 pos;
  Vec2 vel;
  unsigned id = 0;
  bool dead = true;
};

// Free-list of available particle slots
class ParticleSlots {
public:
  ParticleSlots() = default;

  void push(int idx) { free_.push_back(idx); }
  bool empty() const { return free_.empty(); }
  int size() const { return static_cast<int>(free_.size()); }

  // Pop one free slot.  Returns -1 if the list is empty.
  int pop() {
    if (free_.empty())
      return -1;
    int idx = free_.back();
    free_.pop_back();
    return idx;
  }

  // Legacy 2-arg interface used by AdvectParticles / Seeding.
  void AddParticleSlot(int idx, int /*unused*/) { push(idx); }
  bool Empty() const { return empty(); }
  bool PopParticleSlot(int &idx, int &unused) {
    unused = 0;
    idx = pop();
    return idx >= 0;
  }

private:
  std::vector<int> free_;
};

// Particle array
class Particles {
public:
  int nx, ny;
  varType dx, dy;
  int ppcx, ppcy;

  // Total allocated slots.  Always >= nx*ny*ppcx*ppcy.
  int capacity;

  // Kept for backward compatibility with loops that use px/py as bounds.
  // px = capacity, py = 1  →  Index(i,j) = i  (single flat dimension).
  int px;
  int py;

  std::vector<Particle> A;

  // max_particles: total slot count.  Must be >= nx*ny*ppcx*ppcy.
  // Typically set to 2× or 3× that to provide a reservoir for inflow.
  Particles(int nx, int ny, varType dx, varType dy, int ppcx, int ppcy,
            int max_particles = -1)
      : nx(nx), ny(ny), dx(dx), dy(dy), ppcx(ppcx), ppcy(ppcy),
        capacity(max_particles > 0 ? max_particles : ppcx * nx * ppcy * ny),
        px(capacity), py(1), A(static_cast<std::size_t>(capacity)) {}

  // --- accessors (flat index) ---
  varType GetX(int i) const { return A[i].pos.x; }
  varType GetY(int i) const { return A[i].pos.y; }
  varType GetU(int i) const { return A[i].vel.x; }
  varType GetV(int i) const { return A[i].vel.y; }
  unsigned GetId(int i) const { return A[i].id; }
  bool IsDead(int i) const { return A[i].dead; }

  void SetX(int i, varType v) { A[i].pos.x = v; }
  void SetY(int i, varType v) { A[i].pos.y = v; }
  void SetU(int i, varType v) { A[i].vel.x = v; }
  void SetV(int i, varType v) { A[i].vel.y = v; }
  void SetId(int i, unsigned v) { A[i].id = v; }
  void SetDead(int i, bool v) { A[i].dead = v; }

  // --- legacy 2-arg accessors (ip unused, jp = flat index) ---
  varType GetX(int /*ip*/, int jp) const { return GetX(jp); }
  varType GetY(int /*ip*/, int jp) const { return GetY(jp); }
  varType GetU(int /*ip*/, int jp) const { return GetU(jp); }
  varType GetV(int /*ip*/, int jp) const { return GetV(jp); }
  unsigned GetId(int /*ip*/, int jp) const { return GetId(jp); }
  bool IsDead(int /*ip*/, int jp) const { return IsDead(jp); }

  void SetX(int /*ip*/, int jp, varType v) { SetX(jp, v); }
  void SetY(int /*ip*/, int jp, varType v) { SetY(jp, v); }
  void SetU(int /*ip*/, int jp, varType v) { SetU(jp, v); }
  void SetV(int /*ip*/, int jp, varType v) { SetV(jp, v); }
  void SetId(int /*ip*/, int jp, varType v) {
    SetId(jp, static_cast<unsigned>(v));
  }
  void SetDead(int /*ip*/, int jp, bool v) { SetDead(jp, v); }

  // Activate slot idx with the given state.
  void DropOneParticle(int idx, varType x, varType y, varType u, varType v,
                       unsigned id);

  // Legacy 3-arg signature (ip ignored, jp = flat index).
  void DropOneParticle(int /*ip*/, int jp, varType x, varType y, varType u,
                       varType v, int id) {
    DropOneParticle(jp, x, y, u, v, static_cast<unsigned>(id));
  }

  // Place the first nx*ny*ppcx*ppcy slots on a regular sub-cell grid.
  // All remaining slots stay dead; the caller should push them onto the
  // free-list (see PIC::InitFreeSlots).
  void InitParticleGrid();
};
