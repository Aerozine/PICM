#include "Fields2D.hpp"

Grid2D Field2D::Div(){

  Grid2D div(nx - 1, ny - 1);
  
  for(size_t j = 0; j < ny - 1; j++){
    for(size_t i = 0; i < nx - 1; i++){
      varType dudx = (u.Get(i, j) - u.Get(i, j + 1)) / nx;
      varType dvdy = (v.Get(i, j) - v.Get(i + 1, j)) / ny;

      div.Set(i, j, dudx + dvdy);
    }
  } 
  return div;
}
