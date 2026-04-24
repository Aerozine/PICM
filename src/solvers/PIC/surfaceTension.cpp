#include "PIC.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

// Steps d'implémentation

// add phi field
// add params.surfaceTension && params.gamma && params.particleRadius
// modify neighbourSum (if voisin is Air -> p += $\gamma \kappa$)
// modify Step() to compute $\phi$

/*
void PIC::UpdatePhiFromParticles() const {
  varType particleRadius = params.particleRadius;
  //@todo may use inf ?
  const varType far = static_cast<varType>(1e6);

  OMP_PRAGMA(omp parallel for collapse(2))
  for (int i = 0; i < fields->phi.nx; ++i) {
    for (int j = 0; j < fields->phi.ny; ++j) {

      // fields->phi has the same indexing as normVelocityCenter / div.
      // Their centers are at:
      // x = (i + 0.5) dx
      // y = (j + 0.5) dy

      const varType x = (static_cast<varType>(i) + static_cast<varType>(0.5)) *
dx; const varType y = (static_cast<varType>(j) + static_cast<varType>(0.5)) *
dy;

      varType minDist = far;

      for (int p = 0; p < particles->size(); ++p) {
        const varType px = particles->GetX(p);
        const varType py = particles->GetY(p);

        const varType rx = x - px;
        const varType ry = y - py;

        const varType dist = std::sqrt(rx * rx + ry * ry) - particleRadius;
        minDist = std::min(minDist, dist);
      }

      // force the sign to agree with your cell labels.
      // phi < 0 : fluid
      // phi > 0 : air / solid

      varType phi = minDist;

      const labeltype label = fields->Label(i, j);

      if (IS_FLUID(label)) {
        phi = -std::abs(phi);
      } else {
        phi = std::abs(phi);
      }

      fields->phi.Set(i, j, phi);
    }
  }
}

// @todo uniform prototype
inline varType SurfaceTensionPressure(int i, int j,const Fields2D &fields,
varType gamma)  { const int im = std::clamp(i - 1, 0, fields.phi.nx - 1); const
int ip = std::clamp(i + 1, 0, fields.phi.nx - 1); const int jm = std::clamp(j -
1, 0, fields.phi.ny - 1); const int jp = std::clamp(j + 1, 0, fields.phi.ny -
1);

  auto normalX = [&](int a, int b) -> varType {
    const int am = std::clamp(a - 1, 0, fields.phi.nx - 1);
    const int ap = std::clamp(a + 1, 0, fields.phi.nx - 1);
    const int bm = std::clamp(b - 1, 0, fields.phi.ny - 1);
    const int bp = std::clamp(b + 1, 0, fields.phi.ny - 1);

    const varType dphidx =
        (fields.phi.Get(ap, b) - fields.phi.Get(am, b)) /
        (static_cast<varType>(2) * fields.dx);

    const varType dphidy =
        (fields.phi.Get(a, bp) - fields.phi.Get(a, bm)) /
        (static_cast<varType>(2) * fields.dy);

    const varType norm =
        std::sqrt(dphidx * dphidx + dphidy * dphidy) +
        static_cast<varType>(1e-12);

    return dphidx / norm;
  };

  auto normalY = [&](int a, int b) -> varType {
    const int am = std::clamp(a - 1, 0, fields.phi.nx - 1);
    const int ap = std::clamp(a + 1, 0, fields.phi.nx - 1);
    const int bm = std::clamp(b - 1, 0, fields.phi.ny - 1);
    const int bp = std::clamp(b + 1, 0, fields.phi.ny - 1);

    const varType dphidx =
        (fields.phi.Get(ap, b) - fields.phi.Get(am, b)) /
        (static_cast<varType>(2) * fields.dx);

    const varType dphidy =
        (fields.phi.Get(a, bp) - fields.phi.Get(a, bm)) /
        (static_cast<varType>(2) * fields.dy);

    const varType norm =
        std::sqrt(dphidx * dphidx + dphidy * dphidy) +
        static_cast<varType>(1e-12);

    return dphidy / norm;
  };

  const varType dnx_dx =
      (normalX(ip, j) - normalX(im, j)) /
      (static_cast<varType>(2) * fields.dx);

  const varType dny_dy =
      (normalY(i, jp) - normalY(i, jm)) /
      (static_cast<varType>(2) * fields.dy);

  const varType kappa = dnx_dx + dny_dy;

  return gamma * kappa;
}
*/