#include "FLIP.hpp"

FLIP::FLIP(Parameters &params)
    : PIC(params), u_old(fields->u.nx, fields->u.ny),
      v_old(fields->v.nx, fields->v.ny) {}

void FLIP::SaveOldVelocities() {
    OMP_PRAGMA(omp parallel for collapse(2))
    for (int j = 0; j < fields->u.ny; ++j)
      for (int i = 0; i < fields->u.nx; ++i)
        u_old.Set(i, j, fields->u.Get(i, j));

    OMP_PRAGMA(omp parallel for collapse(2))
    for (int j = 0; j < fields->v.ny; ++j)
      for (int i = 0; i < fields->v.nx; ++i)
        v_old.Set(i, j, fields->v.Get(i, j));
}

void FLIP::Step() {
  if (params.gravity > 0.0)
    ApplyGravity();
  ProjectParticlesOnGrid();
  SaveOldVelocities(); // save old veloctiy (!= PIC)
  MakeIncompressible(params, *fields);
  ProjectGridOnParticles(); // based on u_old and u_new (!= PIC)
  fields->Div();
  fields->VelocityNormCenterGrid();
  Advect();
  CountAliveParticles();
  UpdateCellState();
  if (params.refill)
    RefillParticles();
}

void FLIP::ProjectGridOnParticles() {
    varType coefPic  = 0.05; 
    varType coefFlip = 0.95; 
    OMP_PRAGMA(omp parallel for)
    for (int idx = 0; idx < particles->size(); ++idx) {
      const varType x = particles->GetX(idx);
      const varType y = particles->GetY(idx);

      const varType up_old = particles->GetU(idx);
      const varType vp_old = particles->GetV(idx);

      const varType u_new_grid = interpolateU(fields->u, x, y);
      const varType v_new_grid = interpolateV(fields->v, x, y);

      const varType u_old_grid = interpolateU(u_old, x, y);
      const varType v_old_grid = interpolateV(v_old, x, y);

      particles->SetU(idx, coefFlip * (up_old + (u_new_grid - u_old_grid))
                   + coefPic  * u_new_grid);
      particles->SetV(idx, coefFlip * (vp_old + (v_new_grid - v_old_grid))
                   + coefPic  * v_new_grid);

      // particles->SetU(idx, up_old + (u_new_grid - u_old_grid));
      // particles->SetV(idx, vp_old + (v_new_grid - v_old_grid));
    }
}