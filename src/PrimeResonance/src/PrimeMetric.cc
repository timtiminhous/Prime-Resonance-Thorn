#include "cctk.h"
    #include "cctk_Arguments.h"
    #include "cctk_Parameters.h"
    #include <cmath>

    extern "C" void PrimeResonance_Init(CCTK_ARGUMENTS) {
        DECLARE_CCTK_ARGUMENTS
        DECLARE_CCTK_PARAMETERS

        int index;
        // FIX: Correctly loop over grid points
        for(int k=0; k<cctk_lsh[2]; k++)
          for(int j=0; j<cctk_lsh[1]; j++)
            for(int i=0; i<cctk_lsh[0]; i++) {
                // FIX: Use cctkGH (Capitalized)
                index = CCTK_GFINDEX3D(cctkGH, i, j, k);
                
                // FIX: Access coordinate arrays with [index]
                CCTK_REAL r2 = x[index]*x[index] + y[index]*y[index] + z[index]*z[index];
                
                phi[index] = initial_amplitude * exp(-r2 / (pulse_width * pulse_width));
                phi_t[index] = 0;
            }
    }

    extern "C" void PrimeResonance_CalcTmunu(CCTK_ARGUMENTS) {
        DECLARE_CCTK_ARGUMENTS
        DECLARE_CCTK_PARAMETERS

        int index;
        for(int k=0; k<cctk_lsh[2]; k++)
          for(int j=0; j<cctk_lsh[1]; j++)
            for(int i=0; i<cctk_lsh[0]; i++) {
                // FIX: Use cctkGH (Capitalized)
                index = CCTK_GFINDEX3D(cctkGH, i, j, k);
                
                
      // 1. Calculate Radius (Safe)
      // FIX: Access arrays with [index]
      CCTK_REAL r = sqrt(x[index]*x[index] + y[index]*y[index] + z[index]*z[index]) + 1.0e-6;
      
      // 2. The Prime Constants
      CCTK_REAL alpha_inv = 137.035999;
      
      // 3. The Waterfall Potential V(phi)
      CCTK_REAL V_prime = sin(phi[index] * alpha_inv) / r + 0.5 * phi[index] * phi[index];
      
      // 4. Kinetic Term
      CCTK_REAL kinetic = 0.5 * phi_t[index] * phi_t[index];
      
      // 5. Compute Stress-Energy
      eTtt[index] = kinetic + V_prime;
      
      // 6. Spatial Stress (Pressure)
      CCTK_REAL pressure = kinetic - V_prime;
      eTxx[index] = pressure;
      eTyy[index] = pressure;
      eTzz[index] = pressure;

            }
    }