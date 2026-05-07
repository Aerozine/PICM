#include "FLIP.hpp"

FLIP::FLIP(Parameters &params)
    : PIC(params), u_old(fields->u.nx, fields->u.ny),
      v_old(fields->v.nx, fields->v.ny) {}

void FLIP::SaveOldVelocities() {
  u_old = fields->u;
  v_old = fields->v;
}

void FLIP::Step() {
  const varType frameDt = dt;
  const int substeps = computeAdvectionSubsteps();
  const varType subDt = frameDt / static_cast<varType>(substeps);

  if (substeps > 1) {
    DBG_PRINTF("%s: CFL substepping %d x dt=%g",
               solverMethodName(params.solver.method).c_str(), substeps,
               static_cast<double>(subDt));
  }

  setTimeStep(subDt);
  for (int substep = 0; substep < substeps; ++substep) {
    ProjectParticlesOnGrid();
    SaveOldVelocities();
    if (params.surfaceTension && !params.particleInteraction)
      LaplacePressure();
    MakeIncompressible(params, *fields);
    ProjectGridOnParticles();
    if (params.surfaceTension) {
      particleInteraction();
      // Surface tension updates particle velocities. Reproject them so the
      // particle advection in this same substep sees the capillary kick.
      ProjectParticlesOnGrid();
    }
    fields->UpdateDivNorm();
    Advect();
    UpdateCellState();
    if (params.refill)
      RefillParticles();
  }
  setTimeStep(frameDt);
}

void FLIP::ProjectGridOnParticles() {
  const varType coefPic = params.coefPic;
  const varType coefFlip = varType(1) - coefPic;

  OMP_PRAGMA(omp parallel for collapse(2) schedule(static))
  for (int cj = 0; cj < ny; ++cj) {
    for (int ci = 0; ci < nx; ++ci) {
      Particles &cell = (*cloud)(ci, cj);
      for (int p = 0; p < cell.size(); ++p) {
        const varType x = cell.GetX(p);
        const varType y = cell.GetY(p);

        const varType up_old = cell.GetU(p);
        const varType vp_old = cell.GetV(p);

        const varType u_new_grid = fields->interpolateU(x, y);
        const varType v_new_grid = fields->interpolateV(x, y);

        const varType u_old_grid = fields->interpolateU(u_old, x, y);
        const varType v_old_grid = fields->interpolateV(v_old, x, y);

        cell.SetU(p, coefPic * u_new_grid +
                         coefFlip * (up_old + (u_new_grid - u_old_grid)));
        cell.SetV(p, coefPic * v_new_grid +
                         coefFlip * (vp_old + (v_new_grid - v_old_grid)) -
                         dt * params.gravity);
      }
    }
  }
}
