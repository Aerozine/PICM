#include "PIC.hpp"

void PIC::particleInteraction() {
    const int R = 1;

    OMP_PRAGMA(omp parallel for collapse(2) schedule(dynamic,1))
    for (int j = 0; j < fields->ny; ++j) {
        for (int i = 0; i < fields->nx; ++i) {
            
            constexpr int interfaceWidth = 3;
             int ci_lo = std::max(0, i - interfaceWidth);
             int ci_hi = std::min(nx - 1, i + interfaceWidth);
             int cj_lo = std::max(0, j - interfaceWidth);
             int cj_hi = std::min(ny - 1, j + interfaceWidth);

            bool airpresent = false;
            for (int cj = cj_lo; cj <= cj_hi; ++cj) {
                    for (int ci = ci_lo; ci <= ci_hi; ++ci) {
            
            const labeltype left   = fields->Label(ci, cj + 1);
            const labeltype right  = fields->Label(ci + 2, cj + 1);
            const labeltype up     = fields->Label(ci + 1, cj);
            const labeltype down   = fields->Label(ci + 1, cj + 2);

            airpresent|=IS_AIR(left) || IS_AIR(right) || IS_AIR(up) || IS_AIR(down);
                    }
            }
			labeltype center = fields->Label(i + 1, j +1);
            if (!airpresent)
                continue;

            if (!IS_FLUID(center))
                continue;

            //const bool isInterface =
            //    IS_AIR(left) || IS_AIR(right) || IS_AIR(up) || IS_AIR(down);

            //if (!isInterface)
            //    continue;
            
            Particles &centerCell = (*cloud)(i , j);
            const int particleCount = centerCell.size();
            constexpr int AAAAA = 3;
             ci_lo = std::max(0, i - AAAAA);
             ci_hi = std::min(nx - 1, i + AAAAA);
             cj_lo = std::max(0, j - AAAAA);
             cj_hi = std::min(ny - 1, j + AAAAA);

            for (int p = 0; p < particleCount; ++p) {
                varType du = 0.0;
                varType dv = 0.0;

                const varType xp = centerCell.GetX(p);
                const varType yp = centerCell.GetY(p);

                for (int cj = cj_lo; cj <= cj_hi; ++cj) {
                    for (int ci = ci_lo; ci <= ci_hi; ++ci) {

                        if (!IS_FLUID(fields->Label(ci + 1, cj + 1)))
                            continue;

                        const Particles &neighborCell = (*cloud)(ci, cj);
                        const int neighborParticleCount = neighborCell.size();

                        for (int q = 0; q < neighborParticleCount; ++q) {
                            if (ci == i && cj == j && q == p)
                                continue;

                            const varType xq = neighborCell.GetX(q);
                            const varType yq = neighborCell.GetY(q);

                            const varType dx = xq - xp;
                            const varType dy = yq - yp;

                            const varType r2 = dx * dx + dy * dy;

                            if (r2 < REAL_EPSILON)
                                continue;

                            const varType r = std::sqrt(r2);

                            const varType interactionRadius =
                                params.particleRadius*10.0 ;

                            if (r > interactionRadius)
                                continue;

                            const varType nx = dx / r;
                            const varType ny = dy / r;

                            const varType strength =
                                params.dt * params.interactionStiffness * 
                                (1.0 - ((r) / (interactionRadius)));

                            du += strength * nx;
                            dv += strength * ny;
                        }
                    }
                }

                centerCell.SetU(p, centerCell.GetU(p) + du);
                centerCell.SetV(p, centerCell.GetV(p) + dv);
            }
        }
    }
}