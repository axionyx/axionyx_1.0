
# File fdm\_3d.f90


[**File List**](files.md) **>** [**FDM**](dir_43b815edcf2a06ee60d8a45cc6c77fb8.md) **>** [**fdm\_3d.f90**](fdm__3d_8f90.md)

[Go to the source code of this file.](fdm__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**deposit\_fdm\_particles**](fdm__3d_8f90.md#function-deposit-fdm-particles) (particles particles, np np, state\_real state\_real, lo\_real lo\_real, hi\_real hi\_real, state\_imag state\_imag, lo\_imag lo\_imag, hi\_imag hi\_imag, plo plo, dx dx, a a) <br> |
|  subroutine | [**deposit\_fdm\_particles\_wkb**](fdm__3d_8f90.md#function-deposit-fdm-particles-wkb) (particles particles, np np, state\_real state\_real, lo\_real lo\_real, hi\_real hi\_real, state\_imag state\_imag, lo\_imag lo\_imag, hi\_imag hi\_imag, plo plo, dx dx, a a) <br> |
|  subroutine | [**fort\_fdm\_fields**](fdm__3d_8f90.md#function-fort-fdm-fields) (state state, state\_l1 state\_l1, state\_l2 state\_l2, state\_l3 state\_l3, state\_h1 state\_h1, state\_h2 state\_h2, state\_h3 state\_h3) <br> |
|  subroutine | [**fort\_set\_a**](fdm__3d_8f90.md#function-fort-set-a) (scalefactor scalefactor) <br> |
|  subroutine | [**fort\_set\_gamma**](fdm__3d_8f90.md#function-fort-set-gamma) (gamma gamma) <br> |
|  subroutine | [**fort\_set\_hbaroverm**](fdm__3d_8f90.md#function-fort-set-hbaroverm) (hbm hbm) <br> |
|  subroutine | [**fort\_set\_meandens**](fdm__3d_8f90.md#function-fort-set-meandens) (dens dens) <br> |
|  subroutine | [**fort\_set\_mtt**](fdm__3d_8f90.md#function-fort-set-mtt) (mtt mtt) <br> |
|  subroutine | [**fort\_set\_sigma**](fdm__3d_8f90.md#function-fort-set-sigma) (sigma sigma) <br> |
|  subroutine | [**fort\_set\_theta**](fdm__3d_8f90.md#function-fort-set-theta) (theta theta) <br> |








## Public Functions Documentation


### <a href="#function-deposit-fdm-particles" id="function-deposit-fdm-particles">function deposit\_fdm\_particles </a>


```cpp
subroutine deposit_fdm_particles (
    particles particles,
    np np,
    state_real state_real,
    lo_real lo_real,
    hi_real hi_real,
    state_imag state_imag,
    lo_imag lo_imag,
    hi_imag hi_imag,
    plo plo,
    dx dx,
    a a
) 
```



### <a href="#function-deposit-fdm-particles-wkb" id="function-deposit-fdm-particles-wkb">function deposit\_fdm\_particles\_wkb </a>


```cpp
subroutine deposit_fdm_particles_wkb (
    particles particles,
    np np,
    state_real state_real,
    lo_real lo_real,
    hi_real hi_real,
    state_imag state_imag,
    lo_imag lo_imag,
    hi_imag hi_imag,
    plo plo,
    dx dx,
    a a
) 
```



### <a href="#function-fort-fdm-fields" id="function-fort-fdm-fields">function fort\_fdm\_fields </a>


```cpp
subroutine fort_fdm_fields (
    state state,
    state_l1 state_l1,
    state_l2 state_l2,
    state_l3 state_l3,
    state_h1 state_h1,
    state_h2 state_h2,
    state_h3 state_h3
) 
```



### <a href="#function-fort-set-a" id="function-fort-set-a">function fort\_set\_a </a>


```cpp
subroutine fort_set_a (
    scalefactor scalefactor
) 
```



### <a href="#function-fort-set-gamma" id="function-fort-set-gamma">function fort\_set\_gamma </a>


```cpp
subroutine fort_set_gamma (
    gamma gamma
) 
```



### <a href="#function-fort-set-hbaroverm" id="function-fort-set-hbaroverm">function fort\_set\_hbaroverm </a>


```cpp
subroutine fort_set_hbaroverm (
    hbm hbm
) 
```



### <a href="#function-fort-set-meandens" id="function-fort-set-meandens">function fort\_set\_meandens </a>


```cpp
subroutine fort_set_meandens (
    dens dens
) 
```



### <a href="#function-fort-set-mtt" id="function-fort-set-mtt">function fort\_set\_mtt </a>


```cpp
subroutine fort_set_mtt (
    mtt mtt
) 
```



### <a href="#function-fort-set-sigma" id="function-fort-set-sigma">function fort\_set\_sigma </a>


```cpp
subroutine fort_set_sigma (
    sigma sigma
) 
```



### <a href="#function-fort-set-theta" id="function-fort-set-theta">function fort\_set\_theta </a>


```cpp
subroutine fort_set_theta (
    theta theta
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/FDM/fdm_3d.f90`