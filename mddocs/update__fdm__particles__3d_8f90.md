
# File update\_fdm\_particles\_3d.f90


[**File List**](files.md) **>** [**FDM**](dir_43b815edcf2a06ee60d8a45cc6c77fb8.md) **>** [**update\_fdm\_particles\_3d.f90**](update__fdm__particles__3d_8f90.md)

[Go to the source code of this file.](update__fdm__particles__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**update\_fdm\_particles**](update__fdm__particles__3d_8f90.md#function-update-fdm-particles) (np np, particles particles, accel accel, accel\_lo accel\_lo, accel\_hi accel\_hi, plo plo, dx dx, dt dt, a\_prev a\_prev, a\_cur a\_cur, do\_move do\_move) <br> |
|  subroutine | [**update\_fdm\_particles\_wkb**](update__fdm__particles__3d_8f90.md#function-update-fdm-particles-wkb) (np np, particles particles, accel accel, accel\_lo accel\_lo, accel\_hi accel\_hi, plo plo, dx dx, dt dt, a\_prev a\_prev, a\_cur a\_cur, do\_move do\_move) <br> |
|  subroutine | [**update\_gaussian\_beams**](update__fdm__particles__3d_8f90.md#function-update-gaussian-beams) (np np, particles particles, accel accel, accel\_lo accel\_lo, accel\_hi accel\_hi, phi phi, phi\_lo phi\_lo, phi\_hi phi\_hi, plo plo, dx dx, dt dt, a\_prev a\_prev, a\_cur a\_cur, do\_move do\_move) <br> |
|  subroutine | [**update\_gaussian\_beams\_wkb**](update__fdm__particles__3d_8f90.md#function-update-gaussian-beams-wkb) (np np, particles particles, accel accel, accel\_lo accel\_lo, accel\_hi accel\_hi, phi phi, phi\_lo phi\_lo, phi\_hi phi\_hi, plo plo, dx dx, dt dt, a\_prev a\_prev, a\_cur a\_cur, do\_move do\_move) <br> |








## Public Functions Documentation


### <a href="#function-update-fdm-particles" id="function-update-fdm-particles">function update\_fdm\_particles </a>


```cpp
subroutine update_fdm_particles (
    np np,
    particles particles,
    accel accel,
    accel_lo accel_lo,
    accel_hi accel_hi,
    plo plo,
    dx dx,
    dt dt,
    a_prev a_prev,
    a_cur a_cur,
    do_move do_move
) 
```



### <a href="#function-update-fdm-particles-wkb" id="function-update-fdm-particles-wkb">function update\_fdm\_particles\_wkb </a>


```cpp
subroutine update_fdm_particles_wkb (
    np np,
    particles particles,
    accel accel,
    accel_lo accel_lo,
    accel_hi accel_hi,
    plo plo,
    dx dx,
    dt dt,
    a_prev a_prev,
    a_cur a_cur,
    do_move do_move
) 
```



### <a href="#function-update-gaussian-beams" id="function-update-gaussian-beams">function update\_gaussian\_beams </a>


```cpp
subroutine update_gaussian_beams (
    np np,
    particles particles,
    accel accel,
    accel_lo accel_lo,
    accel_hi accel_hi,
    phi phi,
    phi_lo phi_lo,
    phi_hi phi_hi,
    plo plo,
    dx dx,
    dt dt,
    a_prev a_prev,
    a_cur a_cur,
    do_move do_move
) 
```



### <a href="#function-update-gaussian-beams-wkb" id="function-update-gaussian-beams-wkb">function update\_gaussian\_beams\_wkb </a>


```cpp
subroutine update_gaussian_beams_wkb (
    np np,
    particles particles,
    accel accel,
    accel_lo accel_lo,
    accel_hi accel_hi,
    phi phi,
    phi_lo phi_lo,
    phi_hi phi_hi,
    plo plo,
    dx dx,
    dt dt,
    a_prev a_prev,
    a_cur a_cur,
    do_move do_move
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/FDM/update_fdm_particles_3d.f90`