#include "FLIP.hpp"

FLIP::FLIP(Parameters &params)
    : PIC(params), u_old(fields->u.nx, fields->u.ny),
      v_old(fields->v.nx, fields->v.ny) {}

void FLIP::SaveOldVelocities() {
  u_old = fields->u;
  v_old = fields->v;
}

void FLIP::Step() {
  ProjectParticlesOnGrid();
  SaveOldVelocities();
  MakeIncompressible(params, *fields);
  ProjectGridOnParticles();
  fields->UpdateDivNorm();
  Advect();
  UpdateCellState();
  if (params.refill)
    RefillParticles();
}

void FLIP::ProjectGridOnParticles() {
  const varType coefPic = params.coefPic;
  const varType coefFlip = varType(1) - coefPic;

  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int ci = 0; ci < nx; ++ci) {
    for (int cj = 0; cj < ny; ++cj) {
      Particles &cell = (*cloud)(ci, cj);
      for (int p = 0; p < cell.size(); ++p) {
        const varType x = cell.GetX(p);
        const varType y = cell.GetY(p);

        const varType up_old = cell.GetU(p);
        const varType vp_old = cell.GetV(p);

        const varType u_new_grid = fields->u.interpolate<0>(x, y, dx, dy);
        const varType v_new_grid = fields->v.interpolate<1>(x, y, dx, dy);

        const varType u_old_grid = u_old.interpolate<0>(x, y, dx, dy);
        const varType v_old_grid = v_old.interpolate<1>(x, y, dx, dy);

        cell.SetU(p, coefPic * u_new_grid +
                         coefFlip * (up_old + (u_new_grid - u_old_grid)));
        cell.SetV(p, coefPic * v_new_grid +
                         coefFlip * (vp_old + (v_new_grid - v_old_grid)) -
                         dt * params.gravity);
      }
    }
  }
}
