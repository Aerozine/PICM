#pragma once
#include "Precision.hpp"
#include <vector>

/**
 * @file Particles.hpp
 * @brief Storage of the properties of every particule used in the PIC solver.
 */

struct Vec2 {
    varType x = 0.0;
    varType y = 0.0;
};

struct Particle {
    Vec2 x; 
    Vec2 v; 
    int id;
};

class Particles {
public:
  int nx; 
  int ny; 
  varType dx; 
  varType dy;
  
  int ppcx; 
  int ppcy; 

  int px = ppcx * nx;   
  int py = ppcy * ny;   

  std::vector<Particle> A;

  Particles(int nx, int ny, varType dx, varType dy, int ppcx, int ppcy)
      : nx(nx), ny(ny), dx(dx), dy(dy), 
        ppcx(ppcx), ppcy(ppcy), A(px * py) {}

  varType GetX(int i, int j) const {return A[Index(i, j)].x.x;}
  varType GetY(int i, int j) const {return A[Index(i, j)].x.y;}
  varType GetU(int i, int j) const {return A[Index(i, j)].v.x;}
  varType GetV(int i, int j) const {return A[Index(i, j)].v.y;}
  varType GetId(int i, int j) const {return A[Index(i, j)].id;}
  
  void SetX(int i, int j, varType val) {A[Index(i, j)].x.x = val;}
  void SetY(int i, int j, varType val) {A[Index(i, j)].x.y = val;}
  void SetU(int i, int j, varType val) {A[Index(i, j)].v.x = val;}
  void SetV(int i, int j, varType val) {A[Index(i, j)].v.y = val;}
  void SetId(int i, int j, varType val) {A[Index(i, j)].id = val;}

  void InitParticleGrid();

private:

  int Index(int i, int j) const {return py * i + j;}


};


