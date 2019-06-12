
# File fdm\_F.H


[**File List**](files.md) **>** [**FDM**](dir_43b815edcf2a06ee60d8a45cc6c77fb8.md) **>** [**fdm\_F.H**](fdm__F_8H.md)

[Go to the source code of this file.](fdm__F_8H_source.md)



* `#include <AMReX_BLFort.H>`















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**deposit\_fdm\_particles**](fdm__F_8H.md#function-deposit-fdm-particles) (const void \* particles, const long \* np, const amrex::Real \* state\_real, const int \* lo\_real, const int \* hi\_real, const amrex::Real \* state\_imag, const int \* lo\_imag, const int \* hi\_imag, const amrex::Real \* prob\_lo, const amrex::Real \* dx, const amrex::Real & a) <br> |
|  void | [**deposit\_fdm\_particles\_wkb**](fdm__F_8H.md#function-deposit-fdm-particles-wkb) (const void \* particles, const long \* np, const amrex::Real \* state\_real, const int \* lo\_real, const int \* hi\_real, const amrex::Real \* state\_imag, const int \* lo\_imag, const int \* hi\_imag, const amrex::Real \* prob\_lo, const amrex::Real \* dx, const amrex::Real & a) <br> |
|  void | [**fort\_set\_a**](fdm__F_8H.md#function-fort-set-a) (const amrex::Real & a) <br> |
|  void | [**fort\_set\_gamma**](fdm__F_8H.md#function-fort-set-gamma) (const amrex::Real & gamma\_fdm) <br> |
|  void | [**fort\_set\_hbaroverm**](fdm__F_8H.md#function-fort-set-hbaroverm) (const amrex::Real & hbaroverm) <br> |
|  void | [**fort\_set\_meandens**](fdm__F_8H.md#function-fort-set-meandens) (const amrex::Real & meandens) <br> |
|  void | [**fort\_set\_mtt**](fdm__F_8H.md#function-fort-set-mtt) (const amrex::Real & m\_tt) <br> |
|  void | [**fort\_set\_sigma**](fdm__F_8H.md#function-fort-set-sigma) (const amrex::Real & sigma\_fdm) <br> |
|  void | [**fort\_set\_theta**](fdm__F_8H.md#function-fort-set-theta) (const amrex::Real & theta\_fdm) <br> |
|  void | [**initpart**](fdm__F_8H.md#function-initpart) (const int \* level, const double \* time, const int \* lo, const int \* hi, const int \* nd, BL\_FORT\_FAB\_ARG(dat), const double \* delta, const double \* xlo, const double \* xhi) <br> |
|  void | [**update\_fdm\_particles**](fdm__F_8H.md#function-update-fdm-particles) (const int \* np, void \* particles, const amrex::Real \* accel, const int \* accel\_lo, const int \* accel\_hi, const amrex::Real \* prob\_lo, const amrex::Real \* dx, const amrex::Real & dt, const amrex::Real & a\_prev, const amrex::Real & a\_cur, const int \* do\_move) <br> |
|  void | [**update\_fdm\_particles\_wkb**](fdm__F_8H.md#function-update-fdm-particles-wkb) (const int \* np, void \* particles, const amrex::Real \* accel, const int \* accel\_lo, const int \* accel\_hi, const amrex::Real \* prob\_lo, const amrex::Real \* dx, const amrex::Real & dt, const amrex::Real & a\_prev, const amrex::Real & a\_cur, const int \* do\_move) <br> |
|  void | [**update\_gaussian\_beams**](fdm__F_8H.md#function-update-gaussian-beams) (const int \* np, void \* particles, const amrex::Real \* accel, const int \* accel\_lo, const int \* accel\_hi, const amrex::Real \* phi, const int \* phi\_lo, const int \* phi\_hi, const amrex::Real \* prob\_lo, const amrex::Real \* dx, const amrex::Real & dt, const amrex::Real & a\_prev, const amrex::Real & a\_cur, const int \* do\_move) <br> |
|  void | [**update\_gaussian\_beams\_wkb**](fdm__F_8H.md#function-update-gaussian-beams-wkb) (const int \* np, void \* particles, const amrex::Real \* accel, const int \* accel\_lo, const int \* accel\_hi, const amrex::Real \* phi, const int \* phi\_lo, const int \* phi\_hi, const amrex::Real \* prob\_lo, const amrex::Real \* dx, const amrex::Real & dt, const amrex::Real & a\_prev, const amrex::Real & a\_cur, const int \* do\_move) <br> |








## Public Functions Documentation


### <a href="#function-deposit-fdm-particles" id="function-deposit-fdm-particles">function deposit\_fdm\_particles </a>


```cpp
void deposit_fdm_particles (
    const void * particles,
    const long * np,
    const amrex::Real * state_real,
    const int * lo_real,
    const int * hi_real,
    const amrex::Real * state_imag,
    const int * lo_imag,
    const int * hi_imag,
    const amrex::Real * prob_lo,
    const amrex::Real * dx,
    const amrex::Real & a
) 
```



### <a href="#function-deposit-fdm-particles-wkb" id="function-deposit-fdm-particles-wkb">function deposit\_fdm\_particles\_wkb </a>


```cpp
void deposit_fdm_particles_wkb (
    const void * particles,
    const long * np,
    const amrex::Real * state_real,
    const int * lo_real,
    const int * hi_real,
    const amrex::Real * state_imag,
    const int * lo_imag,
    const int * hi_imag,
    const amrex::Real * prob_lo,
    const amrex::Real * dx,
    const amrex::Real & a
) 
```



### <a href="#function-fort-set-a" id="function-fort-set-a">function fort\_set\_a </a>


```cpp
void fort_set_a (
    const amrex::Real & a
) 
```



### <a href="#function-fort-set-gamma" id="function-fort-set-gamma">function fort\_set\_gamma </a>


```cpp
void fort_set_gamma (
    const amrex::Real & gamma_fdm
) 
```



### <a href="#function-fort-set-hbaroverm" id="function-fort-set-hbaroverm">function fort\_set\_hbaroverm </a>


```cpp
void fort_set_hbaroverm (
    const amrex::Real & hbaroverm
) 
```



### <a href="#function-fort-set-meandens" id="function-fort-set-meandens">function fort\_set\_meandens </a>


```cpp
void fort_set_meandens (
    const amrex::Real & meandens
) 
```



### <a href="#function-fort-set-mtt" id="function-fort-set-mtt">function fort\_set\_mtt </a>


```cpp
void fort_set_mtt (
    const amrex::Real & m_tt
) 
```



### <a href="#function-fort-set-sigma" id="function-fort-set-sigma">function fort\_set\_sigma </a>


```cpp
void fort_set_sigma (
    const amrex::Real & sigma_fdm
) 
```



### <a href="#function-fort-set-theta" id="function-fort-set-theta">function fort\_set\_theta </a>


```cpp
void fort_set_theta (
    const amrex::Real & theta_fdm
) 
```



### <a href="#function-initpart" id="function-initpart">function initpart </a>


```cpp
void initpart (
    const int * level,
    const double * time,
    const int * lo,
    const int * hi,
    const int * nd,
    BL_FORT_FAB_ARG(dat),
    const double * delta,
    const double * xlo,
    const double * xhi
) 
```



### <a href="#function-update-fdm-particles" id="function-update-fdm-particles">function update\_fdm\_particles </a>


```cpp
void update_fdm_particles (
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



### <a href="#function-update-fdm-particles-wkb" id="function-update-fdm-particles-wkb">function update\_fdm\_particles\_wkb </a>


```cpp
void update_fdm_particles_wkb (
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



### <a href="#function-update-gaussian-beams" id="function-update-gaussian-beams">function update\_gaussian\_beams </a>


```cpp
void update_gaussian_beams (
    const int * np,
    void * particles,
    const amrex::Real * accel,
    const int * accel_lo,
    const int * accel_hi,
    const amrex::Real * phi,
    const int * phi_lo,
    const int * phi_hi,
    const amrex::Real * prob_lo,
    const amrex::Real * dx,
    const amrex::Real & dt,
    const amrex::Real & a_prev,
    const amrex::Real & a_cur,
    const int * do_move
) 
```



### <a href="#function-update-gaussian-beams-wkb" id="function-update-gaussian-beams-wkb">function update\_gaussian\_beams\_wkb </a>


```cpp
void update_gaussian_beams_wkb (
    const int * np,
    void * particles,
    const amrex::Real * accel,
    const int * accel_lo,
    const int * accel_hi,
    const amrex::Real * phi,
    const int * phi_lo,
    const int * phi_hi,
    const amrex::Real * prob_lo,
    const amrex::Real * dx,
    const amrex::Real & dt,
    const amrex::Real & a_prev,
    const amrex::Real & a_cur,
    const int * do_move
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/FDM/fdm_F.H`