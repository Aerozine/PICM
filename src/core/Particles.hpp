#pragma once
#include "Fields.hpp"
#include "Precision.hpp"
#include <vector>

struct Vec2 {
  varType x = 0.0;
  varType y = 0.0;
};

struct Particle {
  Vec2 pos;
  Vec2 vel;
  Vec2 cu;
  Vec2 cv;
  unsigned id = 0;
};

class Particles {
public:
  int nx, ny, ppcx, ppcy;
  varType dx, dy;
  std::vector<Particle> A;

  Particles(int nx, int ny, varType dx, varType dy, int ppcx, int ppcy)
      : nx(nx), ny(ny), ppcx(ppcx), ppcy(ppcy), dx(dx), dy(dy) {
    // A.reserve(nx * ny * ppcx * ppcy * 3);
  }
  int size() const { return static_cast<int>(A.size()); }
  
  void Add(varType x, varType y, varType u, varType v, unsigned id,
          varType cuX = 0, varType cuY = 0, varType cvX = 0, varType cvY = 0) {
      Particle p;
      p.pos.x = x;
      p.pos.y = y;
      p.vel.x = u;
      p.vel.y = v;
      p.cu.x = cuX;
      p.cu.y = cuY;
      p.cv.x = cvX;
      p.cv.y = cvY;
      p.id = id;
      A.push_back(p);
  }

  void Remove(int i) {
    A[i] = A.back();
    A.pop_back();
  }

  varType GetX(int i) const { return A[i].pos.x; }
  varType GetY(int i) const { return A[i].pos.y; }
  varType GetU(int i) const { return A[i].vel.x; }
  varType GetV(int i) const { return A[i].vel.y; }
  varType GetCuX(int i) const { return A[i].cu.x; }
  varType GetCuY(int i) const { return A[i].cu.y; }
  varType GetCvX(int i) const { return A[i].cv.x; }
  varType GetCvY(int i) const { return A[i].cv.y; }

  void SetX(int i, varType v) { A[i].pos.x = v; }
  void SetY(int i, varType v) { A[i].pos.y = v; }
  void SetU(int i, varType v) { A[i].vel.x = v; }
  void SetV(int i, varType v) { A[i].vel.y = v; }
  void GetCuX(int i, varType v) { A[i].cu.x = v; } 
  void GetCuY(int i, varType v) { A[i].cu.y = v; }
  void GetCvX(int i, varType v) { A[i].cv.x = v; }
  void GetCvY(int i, varType v) { A[i].cv.y = v; }

  void InitParticleGrid(Fields2D &fields);
};