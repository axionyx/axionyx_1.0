
# File agn\_F.H


[**File List**](files.md) **>** [**AGN**](dir_ae7083928535d9dc761b73e4a2ad022f.md) **>** [**agn\_F.H**](agn__F_8H.md)

[Go to the source code of this file.](agn__F_8H_source.md)



* `#include <AMReX_BLFort.H>`















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**agn\_accrete\_mass**](agn__F_8H.md#function-agn-accrete-mass) (const int \* np, void \* particles, const amrex::Real \* state, const amrex::Real \* density\_lost, const int \* s\_lo, const int \* s\_hi, const amrex::Real \* dt, const amrex::Real \* dx) <br> |
|  void | [**agn\_merge\_particles**](agn__F_8H.md#function-agn-merge-particles) (const int \* np, void \* particles, const int \* ng, void \* ghosts, const amrex::Real \* dx) <br> |
|  void | [**agn\_particle\_velocity**](agn__F_8H.md#function-agn-particle-velocity) (const int \* np, void \* particles, const amrex::Real \* state\_old, const int \* sold\_lo, const int \* sold\_hi, const amrex::Real \* state\_new, const int \* snew\_lo, const int \* snew\_hi, const amrex::Real \* dx, const int \* add\_energy) <br> |
|  void | [**agn\_release\_energy**](agn__F_8H.md#function-agn-release-energy) (const int \* np, void \* particles, const amrex::Real \* state, const int \* s\_lo, const int \* s\_hi, const amrex::Real \* diag\_eos, const int \* d\_lo, const int \* d\_hi, const amrex::Real \* a, const amrex::Real \* dx) <br> |
|  void | [**init\_uniform01\_rng**](agn__F_8H.md#function-init-uniform01-rng) () <br> |
|  void | [**nyx\_compute\_overlap**](agn__F_8H.md#function-nyx-compute-overlap) (const int \* np, void \* particles, const int \* ng, void \* ghosts, const amrex::Real \* dx) <br> |
|  void | [**update\_agn\_particles**](agn__F_8H.md#function-update-agn-particles) (const int \* np, void \* particles, const amrex::Real \* accel, const int \* accel\_lo, const int \* accel\_hi, const amrex::Real \* prob\_lo, const amrex::Real \* dx, const amrex::Real & dt, const amrex::Real & a\_prev, const amrex::Real & a\_cur, const int \* do\_move) <br> |








## Public Functions Documentation


### <a href="#function-agn-accrete-mass" id="function-agn-accrete-mass">function agn\_accrete\_mass </a>


```cpp
void agn_accrete_mass (
    const int * np,
    void * particles,
    const amrex::Real * state,
    const amrex::Real * density_lost,
    const int * s_lo,
    const int * s_hi,
    const amrex::Real * dt,
    const amrex::Real * dx
) 
```



### <a href="#function-agn-merge-particles" id="function-agn-merge-particles">function agn\_merge\_particles </a>


```cpp
void agn_merge_particles (
    const int * np,
    void * particles,
    const int * ng,
    void * ghosts,
    const amrex::Real * dx
) 
```



### <a href="#function-agn-particle-velocity" id="function-agn-particle-velocity">function agn\_particle\_velocity </a>


```cpp
void agn_particle_velocity (
    const int * np,
    void * particles,
    const amrex::Real * state_old,
    const int * sold_lo,
    const int * sold_hi,
    const amrex::Real * state_new,
    const int * snew_lo,
    const int * snew_hi,
    const amrex::Real * dx,
    const int * add_energy
) 
```



### <a href="#function-agn-release-energy" id="function-agn-release-energy">function agn\_release\_energy </a>


```cpp
void agn_release_energy (
    const int * np,
    void * particles,
    const amrex::Real * state,
    const int * s_lo,
    const int * s_hi,
    const amrex::Real * diag_eos,
    const int * d_lo,
    const int * d_hi,
    const amrex::Real * a,
    const amrex::Real * dx
) 
```



### <a href="#function-init-uniform01-rng" id="function-init-uniform01-rng">function init\_uniform01\_rng </a>


```cpp
void init_uniform01_rng () 
```



### <a href="#function-nyx-compute-overlap" id="function-nyx-compute-overlap">function nyx\_compute\_overlap </a>


```cpp
void nyx_compute_overlap (
    const int * np,
    void * particles,
    const int * ng,
    void * ghosts,
    const amrex::Real * dx
) 
```



### <a href="#function-update-agn-particles" id="function-update-agn-particles">function update\_agn\_particles </a>


```cpp
void update_agn_particles (
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
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/AGN/agn_F.H`