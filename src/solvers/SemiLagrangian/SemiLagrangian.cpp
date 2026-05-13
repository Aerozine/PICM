#include "SemiLagrangian.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

namespace {

[[nodiscard]] inline bool isSolidSmokeCell(const Fields2D &fields, int i,
                                           int j) noexcept {
  return IS_SOLID(fields.Label(i + 1, j + 1));
}

[[nodiscard]] varType sampleSmokeFreeSpace(const Grid2D &smoke,
                                           const Fields2D &fields, varType x,
                                           varType y, varType dx, varType dy) {
  struct WeightedSample {
    int i;
    int j;
    varType weight;
  };

  const varType iClamped =
      std::clamp(x / dx - REAL_LITERAL(0.5), REAL_LITERAL(0.0),
                 static_cast<varType>(smoke.nx - 1));
  const varType jClamped =
      std::clamp(y / dy - REAL_LITERAL(0.5), REAL_LITERAL(0.0),
                 static_cast<varType>(smoke.ny - 1));

  const int i0 = static_cast<int>(std::floor(iClamped));
  const int j0 = static_cast<int>(std::floor(jClamped));
  const int i1 = std::min(i0 + 1, smoke.nx - 1);
  const int j1 = std::min(j0 + 1, smoke.ny - 1);

  const varType fx = iClamped - static_cast<varType>(i0);
  const varType fy = jClamped - static_cast<varType>(j0);

  const std::array<WeightedSample, 4> samples = {{
      {i0, j0, (REAL_LITERAL(1.0) - fx) * (REAL_LITERAL(1.0) - fy)},
      {i1, j0, fx * (REAL_LITERAL(1.0) - fy)},
      {i0, j1, (REAL_LITERAL(1.0) - fx) * fy},
      {i1, j1, fx * fy},
  }};

  varType weightedSmoke = REAL_LITERAL(0.0);
  varType weightSum = REAL_LITERAL(0.0);
  for (const auto &sample : samples) {
    if (sample.weight <= REAL_LITERAL(0.0) ||
        isSolidSmokeCell(fields, sample.i, sample.j))
      continue;
    weightedSmoke += sample.weight * smoke.Get(sample.i, sample.j);
    weightSum += sample.weight;
  }

  if (weightSum <= REAL_EPSILON)
    return REAL_LITERAL(0.0);
  return weightedSmoke / weightSum;
}

} // namespace

SemiLagrangian::SemiLagrangian(Parameters &params)
    : Solver(params), smokeMap(std::make_unique<Grid2D>(params.nx, params.ny)) {
#ifndef NDEBUG
  std::cout << "Grid dimensions:\n"
            << "  p  (nx,   ny  ): " << fields->p.nx << " x " << fields->p.ny
            << '\n'
            << "  u  (nx+1, ny  ): " << fields->u.nx << " x " << fields->u.ny
            << '\n'
            << "  v  (nx,   ny+1): " << fields->v.nx << " x " << fields->v.ny
            << '\n';
#endif
  params.applyToFields(*fields);
  params.applySmoke(*smokeMap, *fields);
  if (params.write_smoke)
    smokeWriter = std::make_unique<OutputWriter>(params.folder, "smoke");
#ifndef NDEBUG
  std::cout << "SemiLagrangian initialised: " << nx << " x " << ny << " grid, "
            << params.nt << " time steps.\n";
#endif
}

void SemiLagrangian::AdvectSmoke() const {
  Grid2D smokeNew(smokeMap->nx, smokeMap->ny);

    OMP_PRAGMA(omp parallel for collapse(2))
    for (int j = 0; j < smokeMap->ny; ++j) {
      for (int i = 0; i < smokeMap->nx; ++i) {
        if (isSolidSmokeCell(*fields, i, j)) {
          smokeNew.Set(i, j, REAL_LITERAL(0.0));
          continue;
        }

        // if (fields->Label(i, j) & Fields2D::BC_S) {
        if (IS_BC_S(fields->Label(i + 1, j + 1))) {
          assert(std::isfinite(smokeMap->Get(i, j)));
          smokeNew.Set(i, j, smokeMap->Get(i, j));
          continue;
        }

        // Position physique du centre de la cellule (i, j)
        const varType x0 = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * dx;
        const varType y0 = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * dy;

        // RK2 backward trace
        varType u0, v0;
        getVelocity(x0, y0, u0, v0);
        const varType xMid = x0 - REAL_LITERAL(0.5) * dt * u0;
        const varType yMid = y0 - REAL_LITERAL(0.5) * dt * v0;

        varType uMid, vMid;
        getVelocity(xMid, yMid, uMid, vMid);
        varType xDep = x0 - dt * uMid;
        varType yDep = y0 - dt * vMid;

        xDep = std::clamp(xDep, REAL_LITERAL(0.5) * dx,
                          (static_cast<varType>(nx) - REAL_LITERAL(0.5)) * dx);
        yDep = std::clamp(yDep, REAL_LITERAL(0.5) * dy,
                          (static_cast<varType>(ny) - REAL_LITERAL(0.5)) * dy);

        assert(std::isfinite(xDep));
        assert(std::isfinite(yDep));
        smokeNew.Set(
            i, j, sampleSmokeFreeSpace(*smokeMap, *fields, xDep, yDep, dx, dy));
      }
    }
    *smokeMap = std::move(smokeNew);
}

void SemiLagrangian::Step() {
  MakeIncompressible(params, *fields);
  fields->UpdateDivNorm();
  Advect();
  AdvectSmoke();
  DBG_PRINTF("step");
}

void SemiLagrangian::WriteOutput(int step) const {
  if (step % params.sampling_rate != 0)
    return;

  if (params.write_norm_velocity)
    fields->VelocityNormCenterGrid();
  if (params.write_vorticity)
    fields->VorticityCenterGrid();

  Solver::WriteOutput(step);

  if (params.write_smoke && smokeWriter)
    smokeWriter->writeGrid2D(*smokeMap, "smoke",
                             static_cast<double>(step) *
                                 static_cast<double>(dt));
}

void SemiLagrangian::Run() {
  fields->UpdateDivNorm();
  WriteOutput(0);
  RunLoop(std::max(1, params.nt / 20));
}
