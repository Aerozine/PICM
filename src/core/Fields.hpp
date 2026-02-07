#pragma once
#include "Grid2D.hpp"
#include <cstdint> // for uint8_t

class Fields2D {
public:
    typedef float varType;

    enum CellType : uint8_t { FLUID = 0, SOLID = 1 };

    Fields2D(size_t nx, size_t ny, float density)
    : nx(nx), ny(ny), density(density),
      u(nx, ny - 1),
      v(nx - 1, ny),
      p(nx - 1, ny - 1),
      div(nx - 1, ny - 1),
      rot(nx - 1, ny - 1),
      labels(nx * ny, FLUID) {}

    CellType Label(size_t i, size_t j) const {
      return static_cast<CellType>(labels[idx(i,j)]);
    }
    
    void SetLabel(size_t i, size_t j, CellType t) {
      labels[idx(i,j)] = static_cast<uint8_t>(t);
    }

    Grid2D Div();

    bool InBounds(size_t i, size_t j) const {return i < nx && j < ny;}

private:
    float density;
    size_t nx, ny;
  
    varType usolid = 0.0; // first try (fixed borders)
                         // next step -> moving boundaries ? 
    Grid2D u;   
    Grid2D p;  
    Grid2D div;
    Grid2D rot;

    std::vector<uint8_t> labels;
    
    size_t idx(size_t i, size_t j) const {return i + nx * j;}
};

