#include "utilities.hpp"
 
Grid2D div(Grid2D field){

  Grid2D div(field.nx - 1, field.ny - 1);
  
  for(size_t j = 0; j < field.ny - 1; j++){
    for(size_t i = 0; i < field.nx - 1; i++){
      dfxdx = (field.A[i][j] - field.A[i][j + 1])/field.nx;
      dfydy = (field.A[i][j] - field.A[i + 1][j])/field.ny;

      div.A[i][j] = dfxdx + dfydy;
    }
  } 
  return div;
}
