
# File NeutrinoParticles\_F.H

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**NeutrinoParticles\_F.H**](NeutrinoParticles__F_8H.md)

[Go to the documentation of this file.](NeutrinoParticles__F_8H.md) 


````cpp
#ifdef NEUTRINO_PARTICLES
#ifndef _NeutrinoParticles_F_H_
#define _NeutrinoParticles_F_H_

#include <AMReX_BLFort.H>

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    void neutrino_deposit_relativistic_cic(
                                       const amrex_particle_real*, int ns, int np, int nc,
                                       amrex_real* rho, const int* lo, const int* hi,
                                       const amrex_real* plo, const amrex_real* dx,
                                       const amrex_real  csq);

    void neutrino_deposit_particle_dx_relativistic_cic(
                                       const amrex_particle_real*, int ns, int np, int nc,
                                       amrex_real* rho, const int* lo, const int* hi,
                                       const amrex_real* plo, const amrex_real* dx,
                                       const amrex_real* particle_dx,
                                       const amrex_real  csq);

#ifdef __cplusplus
}
#endif

#endif /*NeutrinoParticles_F_H_*/
#endif
````

