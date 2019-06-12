
# File dm\_F.H


[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**dm\_F.H**](dm__F_8H.md)

[Go to the source code of this file.](dm__F_8H_source.md)



* `#include <AMReX_BLFort.H>`















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**update\_dm\_particles**](dm__F_8H.md#function-update-dm-particles) (const int \* np, void \* particles, const amrex::Real \* accel, const int \* accel\_lo, const int \* accel\_hi, const amrex::Real \* prob\_lo, const amrex::Real \* dx, const amrex::Real & dt, const amrex::Real & a\_prev, const amrex::Real & a\_cur, const int \* do\_move) <br> |








## Public Functions Documentation


### <a href="#function-update-dm-particles" id="function-update-dm-particles">function update\_dm\_particles </a>


```cpp
void update_dm_particles (
    const int * np,
    void * particles,
    const amrex::Real * accel,
    const int * accel_lo,
    const int * accel_hi,
    const amrex::Real * prob_lo,
    const amrex::Real * dx,
    const amrex::Real & dt,
    const amrex::Real & a_prev,
    const amrex::Real & a_cur,
    const int * do_move
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/dm_F.H`