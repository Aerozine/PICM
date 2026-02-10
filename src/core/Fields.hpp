#pragma once
#include "Grid2D.hpp"
#include <cstdint> // for uint8_t

class Fields2D {
public:
    typedef float varType;

    enum CellType : uint8_t { FLUID = 0, SOLID = 1 };
    
    size_t nx, ny;
    varType density, dt, dx, dy;
    Grid2D u, v, p, div, rot;   

    Fields2D(size_t nx, size_t ny, varType density,
                  varType dt, varType dx, varType dy)
    : nx(nx), ny(ny), density(density), dt(dt), dx(dx), dy(dy),
      u(nx - 1, ny), // !
      v(nx, ny - 1),
      p(nx - 1, ny - 1),
      div(nx - 1, ny - 1),
      rot(nx - 1, ny - 1),
      labels((nx - 1) * (ny - 1), FLUID) {} 
    
    varType usolid = 0.0; // first try (fixed borders)
                         // next step -> moving boundaries ? 

    CellType Label(size_t i, size_t j) const {
      return static_cast<CellType>(labels[idx(i,j)]);
    }
    
    void SetLabel(size_t i, size_t j, CellType t) {
      labels[idx(i,j)] = static_cast<uint8_t>(t);
    }

    void Div();
    void InitRandomVelocities();
    void InitPotentialGradient(varType amplitude, int kx, int ky);
    //Grid2D VelocityNormCenterGrid();

private:
    std::vector<uint8_t> labels; // not a Grid2D
    size_t idx(size_t i, size_t j) const {return nx * j + i;}
};

