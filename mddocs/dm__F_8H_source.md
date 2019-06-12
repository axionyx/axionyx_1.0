
# File dm\_F.H

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**dm\_F.H**](dm__F_8H.md)

[Go to the documentation of this file.](dm__F_8H.md) 


````cpp
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void update_dm_particles(const int* np, void* particles,
                             const amrex::Real* accel, const int* accel_lo, const int* accel_hi,
                             const amrex::Real* prob_lo, 
                             const amrex::Real* dx, const amrex::Real& dt, 
                             const amrex::Real& a_prev, const amrex::Real& a_cur, const int* do_move);
#ifdef __cplusplus
}
#endif
````

