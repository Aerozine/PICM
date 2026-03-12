#include "PIC.hpp"
#include <iostream>

varType PIC::GetW(){
    int ppcx = particles->ppcx;
    int ppcy = particles->ppcy;

    return ppcx * ppcy; 
}

varType PIC::hat(varType r) {
    if (r >= varType(0) && r <= varType(1))
        return varType(1) - r;
    else if (r >= varType(-1) && r <= varType(0))
        return varType(1) + r;
    else
        return varType(0);
}


void PIC::ProjectOneParticleOnMAC(varType x, varType y, varType up, varType vp)
{
    // projecting on u

    varType xu = x / fields->dx;
    varType yu = y / fields->dy - varType(0.5);

    int i0 = std::floor(xu);
    int j0 = std::floor(yu);

    for (int j = j0; j <= j0 + 1; j++) {
        for (int i = i0; i <= i0 + 1; i++) {
            if (i < 0 || i >= nx + 1 || j < 0 || j >= ny)
                continue;

            varType kx = hat(xu - varType(i));
            varType ky = hat(yu - varType(j));
            varType k  = kx * ky;

            if (k > varType(0))
            {
                fields->u_sum.Set(i, j, k * up);
                fields->u_weight.Set(i, j, k);
            }
        }
    }

    // projecting on v
    varType xv = x / dx - varType(0.5);
    varType yv = y / dy;

    i0 = std::floor(xv);
    j0 = std::floor(yv);

    for (int j = j0; j <= j0 + 1; j++) {
        for (int i = i0; i <= i0 + 1; i++) {
            if (i < 0 || i >= nx || j < 0 || j >= ny + 1)
                continue;

            varType kx = hat(xv - varType(i));
            varType ky = hat(yv - varType(j));
            varType k  = kx * ky;

            if (k > varType(0))
            {
                fields->v_sum.Set(i, j, k * vp);
                fields->v_weight.Set(i, j, k);
            }
        }
    }
}

void PIC::ProjectParticlesOnGrid(std::string kernel)
{
    if (kernel != "hat") {
        std::cout << "invalid kernel for the Particle on Grid projection.\n";
        return;
    }

    // reset all grid values to zero
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx + 1; i++) {
            fields->u_sum.Set(i,j, varType(0)); 
            fields->u_weight.Set(i,j, varType(0));
        }
    }

    for (int j = 1; j < ny + 1; j++) {
        for (int i = 0; i < nx; i++) {
            fields->v_sum.Set(i,j, varType(0)); 
            fields->v_weight.Set(i,j, varType(0));
        }
    }

    int ppcx = particles->ppcx;
    int ppcy = particles->ppcy;

    for (int icell = 0; icell < nx; icell++) {
        for (int jcell = 0; jcell < ny; jcell++) {

            for (int a = 0; a < ppcx; a++) {
                for (int b = 0; b < ppcy; b++) {

                    int ip = icell * ppcx + a;
                    int jp = jcell * ppcy + b;

                    varType x = particles->GetX(ip, jp);
                    varType y = particles->GetY(ip, jp);
                    varType up = particles->GetU(ip, jp);
                    varType vp = particles->GetV(ip, jp);

                    ProjectOneParticleOnMAC(x, y, up, vp);
                }
            }
        }
    }

    varType newVelocity;

    // update u (with normalization)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx + 1; i++) {
            if (fields->u_weight.Get(i,j) > varType(1e-12)) {
                newVelocity =  fields->u_sum.Get(i,j) / fields->u_weight.Get(i,j);
                fields->u.Set(i,j, newVelocity); 
            }
            else
                fields->u.Set(i,j, varType(0));
        }
    }

    // update v (with normalization)
    for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx; i++) {
            if (fields->v_weight.Get(i,j) > varType(1e-12)) {
                newVelocity = fields->v_sum.Get(i,j) / fields->v_weight.Get(i,j);
                fields->v.Set(i,j, newVelocity); 
            }
            else
                fields->v.Set(i,j, varType(0));
        }
    }
}

void PIC::ProjectGridOnParticles(){
    
    int ppcx = particles->ppcx;
    int ppcy = particles->ppcy;

    for (int icell = 0; icell < nx; icell++) {
        for (int jcell = 0; jcell < ny; jcell++) {

            for (int a = 0; a < ppcx; a++) {
                for (int b = 0; b < ppcy; b++) {

                    int ip = icell * ppcx + a;
                    int jp = jcell * ppcy + b;

                    varType x = particles->GetX(ip, jp);
                    varType y = particles->GetY(ip, jp);
                    
                    particles->SetU(ip, jp, interpolateU(x, y));
                    particles->SetV(ip, jp, interpolateV(x, y));
                }
            }
        }
    }
}