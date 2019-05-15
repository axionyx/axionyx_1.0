
# File agn\_3d.f90


[**File List**](files.md) **>** [**AGN**](dir_ae7083928535d9dc761b73e4a2ad022f.md) **>** [**agn\_3d.f90**](agn__3d_8f90.md)

[Go to the source code of this file.](agn__3d_8f90_source.md)












## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**iso\_c\_binding**](namespaceiso__c__binding.md) <br> |






## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**agn\_accrete\_mass**](agn__3d_8f90.md#function-agn-accrete-mass) (np np, particles particles, state state, density\_lost density\_lost, slo slo, shi shi, dt dt, dx dx) <br> |
|  subroutine | [**agn\_merge\_pair**](agn__3d_8f90.md#function-agn-merge-pair) (particle\_stay particle\_stay, particle\_remove particle\_remove) <br> |
|  subroutine | [**agn\_merge\_particles**](agn__3d_8f90.md#function-agn-merge-particles) (np np, particles particles, ng ng, ghosts ghosts, delta\_x delta\_x) <br> |
|  subroutine | [**agn\_particle\_velocity**](agn__3d_8f90.md#function-agn-particle-velocity) (np np, particles particles, state\_old state\_old, sold\_lo sold\_lo, sold\_hi sold\_hi, state\_new state\_new, snew\_lo snew\_lo, snew\_hi snew\_hi, dx dx, add\_energy add\_energy) <br> |
|  subroutine | [**agn\_release\_energy**](agn__3d_8f90.md#function-agn-release-energy) (np np, particles particles, state state, slo slo, shi shi, diag\_eos diag\_eos, dlo dlo, dhi dhi, a a, dx dx) <br> |
|  subroutine | [**get\_length\_frac**](agn__3d_8f90.md#function-get-length-frac) (frac frac, x x, dx dx) <br> |
|  subroutine | [**get\_random\_direction**](agn__3d_8f90.md#function-get-random-direction) (vec vec) <br> |
|  subroutine | [**get\_weights**](agn__3d_8f90.md#function-get-weights) (weight weight, pos pos, dx dx) <br> |
|  subroutine | [**nyx\_compute\_overlap**](agn__3d_8f90.md#function-nyx-compute-overlap) (np np, particles particles, ng ng, ghosts ghosts, delta\_x delta\_x) <br> |








## Public Functions Documentation


### <a href="#function-agn-accrete-mass" id="function-agn-accrete-mass">function agn\_accrete\_mass </a>


```cpp
subroutine agn_accrete_mass (
    np np,
    particles particles,
    state state,
    density_lost density_lost,
    slo slo,
    shi shi,
    dt dt,
    dx dx
) 
```



### <a href="#function-agn-merge-pair" id="function-agn-merge-pair">function agn\_merge\_pair </a>


```cpp
subroutine agn_merge_pair (
    particle_stay particle_stay,
    particle_remove particle_remove
) 
```



### <a href="#function-agn-merge-particles" id="function-agn-merge-particles">function agn\_merge\_particles </a>


```cpp
subroutine agn_merge_particles (
    np np,
    particles particles,
    ng ng,
    ghosts ghosts,
    delta_x delta_x
) 
```



### <a href="#function-agn-particle-velocity" id="function-agn-particle-velocity">function agn\_particle\_velocity </a>


```cpp
subroutine agn_particle_velocity (
    np np,
    particles particles,
    state_old state_old,
    sold_lo sold_lo,
    sold_hi sold_hi,
    state_new state_new,
    snew_lo snew_lo,
    snew_hi snew_hi,
    dx dx,
    add_energy add_energy
) 
```



### <a href="#function-agn-release-energy" id="function-agn-release-energy">function agn\_release\_energy </a>


```cpp
subroutine agn_release_energy (
    np np,
    particles particles,
    state state,
    slo slo,
    shi shi,
    diag_eos diag_eos,
    dlo dlo,
    dhi dhi,
    a a,
    dx dx
) 
```



### <a href="#function-get-length-frac" id="function-get-length-frac">function get\_length\_frac </a>


```cpp
subroutine get_length_frac (
    frac frac,
    x x,
    dx dx
) 
```



### <a href="#function-get-random-direction" id="function-get-random-direction">function get\_random\_direction </a>


```cpp
subroutine get_random_direction (
    vec vec
) 
```



### <a href="#function-get-weights" id="function-get-weights">function get\_weights </a>


```cpp
subroutine get_weights (
    weight weight,
    pos pos,
    dx dx
) 
```



### <a href="#function-nyx-compute-overlap" id="function-nyx-compute-overlap">function nyx\_compute\_overlap </a>


```cpp
subroutine nyx_compute_overlap (
    np np,
    particles particles,
    ng ng,
    ghosts ghosts,
    delta_x delta_x
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/AGN/agn_3d.f90`