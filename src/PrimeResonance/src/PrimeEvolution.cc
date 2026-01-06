#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <cmath>

// Declare MoL function
extern "C" CCTK_INT MoLRegisterEvolved(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex);

extern "C" void PrimeResonance_RegisterVars(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    CCTK_INT ierr = 0;
    
    ierr += MoLRegisterEvolved(CCTK_VarIndex("PrimeResonance::phi"),
                                CCTK_VarIndex("PrimeResonance::phi_rhs"));
    ierr += MoLRegisterEvolved(CCTK_VarIndex("PrimeResonance::phi_t"),
                                CCTK_VarIndex("PrimeResonance::phi_t_rhs"));
    
    if (ierr) CCTK_WARN(0, "Failed to register evolved variables with MoL");
}

extern "C" void PrimeResonance_RHS(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    CCTK_REAL dx = CCTK_DELTA_SPACE(0);
    CCTK_REAL dy = CCTK_DELTA_SPACE(1);
    CCTK_REAL dz = CCTK_DELTA_SPACE(2);
    
    CCTK_REAL idx2 = 1.0 / (dx * dx);
    CCTK_REAL idy2 = 1.0 / (dy * dy);
    CCTK_REAL idz2 = 1.0 / (dz * dz);
    
    CCTK_REAL alpha_inv = 137.035999;
    
    for(int k = 1; k < cctk_lsh[2] - 1; k++)
      for(int j = 1; j < cctk_lsh[1] - 1; j++)
        for(int i = 1; i < cctk_lsh[0] - 1; i++) {
            
            int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
            
            int ip = CCTK_GFINDEX3D(cctkGH, i+1, j, k);
            int im = CCTK_GFINDEX3D(cctkGH, i-1, j, k);
            int jp = CCTK_GFINDEX3D(cctkGH, i, j+1, k);
            int jm = CCTK_GFINDEX3D(cctkGH, i, j-1, k);
            int kp = CCTK_GFINDEX3D(cctkGH, i, j, k+1);
            int km = CCTK_GFINDEX3D(cctkGH, i, j, k-1);
            
            CCTK_REAL laplacian = 
                (phi[ip] - 2.0*phi[index] + phi[im]) * idx2 +
                (phi[jp] - 2.0*phi[index] + phi[jm]) * idy2 +
                (phi[kp] - 2.0*phi[index] + phi[km]) * idz2;
            
            CCTK_REAL r = sqrt(x[index]*x[index] + y[index]*y[index] + z[index]*z[index]) + 1.0e-6;
            
            CCTK_REAL r_safe = sqrt(r*r + 1.0);
            CCTK_REAL dVdphi = alpha_inv * cos(phi[index] * alpha_inv) / r_safe + phi[index];
            
            phi_rhs[index] = phi_t[index];
            phi_t_rhs[index] = laplacian - dVdphi;
        }
    
    for(int k = 0; k < cctk_lsh[2]; k++)
      for(int j = 0; j < cctk_lsh[1]; j++)
        for(int i = 0; i < cctk_lsh[0]; i++) {
            if (i == 0 || i == cctk_lsh[0]-1 ||
                j == 0 || j == cctk_lsh[1]-1 ||
                k == 0 || k == cctk_lsh[2]-1) {
                int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
                phi_rhs[index] = 0.0;
                phi_t_rhs[index] = 0.0;
            }
        }
}
