
# File Nyx\_adv\_FDM\_FD.f90


[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Src\_3d**](dir_723248e6e98dc7cb10ec13b7569a328c.md) **>** [**Nyx\_adv\_FDM\_FD.f90**](Nyx__adv__FDM__FD_8f90.md)

[Go to the source code of this file.](Nyx__adv__FDM__FD_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**fort\_advance\_fdm\_fd**](Nyx__adv__FDM__FD_8f90.md#function-fort-advance-fdm-fd) (time time, lo lo, hi hi, uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1, uin\_h2 uin\_h2, uin\_h3 uin\_h3, uout uout, uout\_l1 uout\_l1, uout\_l2 uout\_l2, uout\_l3 uout\_l3, uout\_h1 uout\_h1, uout\_h2 uout\_h2, uout\_h3 uout\_h3, phi phi, p\_l1 p\_l1, p\_l2 p\_l2, p\_l3 p\_l3, p\_h1 p\_h1, p\_h2 p\_h2, p\_h3 p\_h3, delta delta, prob\_lo prob\_lo, prob\_hi prob\_hi, dt dt, courno courno, a\_old a\_old, a\_new a\_new, verbose verbose) <br>_Solves the Schroedinger-Newton (Gross-Pitaevskii) equation for a self-gravitating Bose-Einstein condensate of an ultralight scalar field using explicit RK4 and directional splitting. (Also, this is an example of how to document Fortran code, in line 300 of Source/Src\_3d/Nyx\_adv\_FDM\_FD.f90)_  |
|  subroutine | [**fort\_divvel**](Nyx__adv__FDM__FD_8f90.md#function-fort-divvel) (uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1) <br> |
|  subroutine | [**fort\_estdtax**](Nyx__adv__FDM__FD_8f90.md#function-fort-estdtax) (uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1, uin\_h2 uin\_h2) <br> |
|  subroutine | [**fort\_initcosmoax**](Nyx__adv__FDM__FD_8f90.md#function-fort-initcosmoax) (init init, i\_l1 i\_l1, i\_l2 i\_l2, i\_l3 i\_l3, i\_h1 i\_h1, i\_h2 i\_h2, i\_h3 i\_h3, state state, s\_l1 s\_l1, s\_l2 s\_l2, s\_l3 s\_l3, s\_h1 s\_h1, s\_h2 s\_h2, s\_h3 s\_h3, divvel divvel, d\_l1 d\_l1, d\_l2 d\_l2, d\_l3 d\_l3, d\_h1 d\_h1, d\_h2 d\_h2, d\_h3 d\_h3, phase phase, p\_l1 p\_l1, p\_l2 p\_l2, p\_l3 p\_l3, p\_h1 p\_h1, p\_h2 p\_h2, p\_h3 p\_h3, didx didx, ivar ivar, dx dx) <br> |
|  subroutine | [**tdma**](Nyx__adv__FDM__FD_8f90.md#function-tdma) (n n, a a, b b, c c, d d, x x) <br> |








## Public Functions Documentation


### <a href="#function-fort-advance-fdm-fd" id="function-fort-advance-fdm-fd">function fort\_advance\_fdm\_fd </a>


```cpp
subroutine fort_advance_fdm_fd (
    time time,
    lo lo,
    hi hi,
    uin uin,
    uin_l1 uin_l1,
    uin_l2 uin_l2,
    uin_l3 uin_l3,
    uin_h1 uin_h1,
    uin_h2 uin_h2,
    uin_h3 uin_h3,
    uout uout,
    uout_l1 uout_l1,
    uout_l2 uout_l2,
    uout_l3 uout_l3,
    uout_h1 uout_h1,
    uout_h2 uout_h2,
    uout_h3 uout_h3,
    phi phi,
    p_l1 p_l1,
    p_l2 p_l2,
    p_l3 p_l3,
    p_h1 p_h1,
    p_h2 p_h2,
    p_h3 p_h3,
    delta delta,
    prob_lo prob_lo,
    prob_hi prob_hi,
    dt dt,
    courno courno,
    a_old a_old,
    a_new a_new,
    verbose verbose
) 
```



### <a href="#function-fort-divvel" id="function-fort-divvel">function fort\_divvel </a>


```cpp
subroutine fort_divvel (
    uin uin,
    uin_l1 uin_l1,
    uin_l2 uin_l2,
    uin_l3 uin_l3,
    uin_h1 uin_h1
) 
```



### <a href="#function-fort-estdtax" id="function-fort-estdtax">function fort\_estdtax </a>


```cpp
subroutine fort_estdtax (
    uin uin,
    uin_l1 uin_l1,
    uin_l2 uin_l2,
    uin_l3 uin_l3,
    uin_h1 uin_h1,
    uin_h2 uin_h2
) 
```



### <a href="#function-fort-initcosmoax" id="function-fort-initcosmoax">function fort\_initcosmoax </a>


```cpp
subroutine fort_initcosmoax (
    init init,
    i_l1 i_l1,
    i_l2 i_l2,
    i_l3 i_l3,
    i_h1 i_h1,
    i_h2 i_h2,
    i_h3 i_h3,
    state state,
    s_l1 s_l1,
    s_l2 s_l2,
    s_l3 s_l3,
    s_h1 s_h1,
    s_h2 s_h2,
    s_h3 s_h3,
    divvel divvel,
    d_l1 d_l1,
    d_l2 d_l2,
    d_l3 d_l3,
    d_h1 d_h1,
    d_h2 d_h2,
    d_h3 d_h3,
    phase phase,
    p_l1 p_l1,
    p_l2 p_l2,
    p_l3 p_l3,
    p_h1 p_h1,
    p_h2 p_h2,
    p_h3 p_h3,
    didx didx,
    ivar ivar,
    dx dx
) 
```



### <a href="#function-tdma" id="function-tdma">function tdma </a>


```cpp
subroutine tdma (
    n n,
    a a,
    b b,
    c c,
    d d,
    x x
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/Src_3d/Nyx_adv_FDM_FD.f90`