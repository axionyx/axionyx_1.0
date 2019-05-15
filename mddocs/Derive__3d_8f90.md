
# File Derive\_3d.f90


[**File List**](files.md) **>** [**DerivedQuantities**](dir_2c61180f16f9dfbd2bd571bcae5f2822.md) **>** [**Derive\_3d.f90**](Derive__3d_8f90.md)

[Go to the source code of this file.](Derive__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**ca\_axangmom\_x**](Derive__3d_8f90.md#function-ca-axangmom-x) (angmom\_x angmom\_x, angmom\_x\_l1 angmom\_x\_l1, angmom\_x\_l2 angmom\_x\_l2) <br> |
|  subroutine | [**ca\_axangmom\_y**](Derive__3d_8f90.md#function-ca-axangmom-y) (angmom\_y angmom\_y, angmom\_y\_l1 angmom\_y\_l1, angmom\_y\_l2 angmom\_y\_l2) <br> |
|  subroutine | [**ca\_axangmom\_z**](Derive__3d_8f90.md#function-ca-axangmom-z) (angmom\_z angmom\_z, angmom\_z\_l1 angmom\_z\_l1, angmom\_z\_l2 angmom\_z\_l2) <br> |
|  subroutine | [**ca\_axekin**](Derive__3d_8f90.md#function-ca-axekin) (ekin ekin, ekin\_l1 ekin\_l1, ekin\_l2 ekin\_l2, ekin\_l3 ekin\_l3, ekin\_h1 ekin\_h1, ekin\_h2 ekin\_h2) <br> |
|  subroutine | [**ca\_axekinrho**](Derive__3d_8f90.md#function-ca-axekinrho) (ekinrho ekinrho, ekinrho\_l1 ekinrho\_l1, ekinrho\_l2 ekinrho\_l2, ekinrho\_l3 ekinrho\_l3) <br> |
|  subroutine | [**ca\_axekinv**](Derive__3d_8f90.md#function-ca-axekinv) (ekinv ekinv, ekinv\_l1 ekinv\_l1, ekinv\_l2 ekinv\_l2, ekinv\_l3 ekinv\_l3, ekinv\_h1 ekinv\_h1) <br> |
|  subroutine | [**ca\_axepot**](Derive__3d_8f90.md#function-ca-axepot) (epot epot, epot\_l1 epot\_l1, epot\_l2 epot\_l2, epot\_l3 epot\_l3, epot\_h1 epot\_h1, epot\_h2 epot\_h2) <br> |
|  subroutine | [**ca\_axphase**](Derive__3d_8f90.md#function-ca-axphase) (phase phase, phase\_l1 phase\_l1, phase\_l2 phase\_l2, phase\_l3 phase\_l3, phase\_h1 phase\_h1) <br> |
|  subroutine | [**ca\_axvel**](Derive__3d_8f90.md#function-ca-axvel) (ekinv ekinv, ekinv\_l1 ekinv\_l1, ekinv\_l2 ekinv\_l2, ekinv\_l3 ekinv\_l3, ekinv\_h1 ekinv\_h1) <br> |
|  subroutine | [**ca\_dererrx**](Derive__3d_8f90.md#function-ca-dererrx) (err err, err\_l1 err\_l1, err\_l2 err\_l2, err\_l3 err\_l3, err\_h1 err\_h1, err\_h2 err\_h2) <br> |
|  subroutine | [**ca\_dererry**](Derive__3d_8f90.md#function-ca-dererry) (err err, err\_l1 err\_l1, err\_l2 err\_l2, err\_l3 err\_l3, err\_h1 err\_h1, err\_h2 err\_h2) <br> |
|  subroutine | [**ca\_dererrz**](Derive__3d_8f90.md#function-ca-dererrz) (err err, err\_l1 err\_l1, err\_l2 err\_l2, err\_l3 err\_l3, err\_h1 err\_h1, err\_h2 err\_h2) <br> |
|  subroutine | [**derdivu**](Derive__3d_8f90.md#function-derdivu) (divu divu, div\_l1 div\_l1, div\_l2 div\_l2, div\_l3 div\_l3, div\_h1 div\_h1, div\_h2 div\_h2, div\_h3 div\_h3) <br> |
|  subroutine | [**dereint1**](Derive__3d_8f90.md#function-dereint1) (e e, e\_l1 e\_l1, e\_l2 e\_l2, e\_l3 e\_l3, e\_h1 e\_h1, e\_h2 e\_h2, e\_h3 e\_h3, ncomp\_e ncomp\_e, u u, u\_l1 u\_l1, u\_l2 u\_l2, u\_l3 u\_l3, u\_h1 u\_h1, u\_h2 u\_h2, u\_h3 u\_h3, ncomp\_u ncomp\_u, lo lo, hi hi, domlo domlo, domhi domhi, dx dx, xlo xlo, time time, dt dt, bc bc, level level, grid\_no grid\_no) <br> |
|  subroutine | [**dereint2**](Derive__3d_8f90.md#function-dereint2) (e e, e\_l1 e\_l1, e\_l2 e\_l2, e\_l3 e\_l3, e\_h1 e\_h1, e\_h2 e\_h2, e\_h3 e\_h3, ncomp\_e ncomp\_e, u u, u\_l1 u\_l1, u\_l2 u\_l2, u\_l3 u\_l3, u\_h1 u\_h1, u\_h2 u\_h2, u\_h3 u\_h3, ncomp\_u ncomp\_u, lo lo, hi hi, domlo domlo, domhi domhi, dx dx, xlo xlo, time time, dt dt, bc bc, level level, grid\_no grid\_no) <br> |
|  subroutine | [**derentropy**](Derive__3d_8f90.md#function-derentropy) (s s, s\_l1 s\_l1, s\_l2 s\_l2, s\_l3 s\_l3, s\_h1 s\_h1, s\_h2 s\_h2, s\_h3 s\_h3, ncomp\_s ncomp\_s, u u, u\_l1 u\_l1, u\_l2 u\_l2, u\_l3 u\_l3, u\_h1 u\_h1, u\_h2 u\_h2, u\_h3 u\_h3, ncomp\_u ncomp\_u) <br> |
|  subroutine | [**derkineng**](Derive__3d_8f90.md#function-derkineng) (kineng kineng, ken\_l1 ken\_l1, ken\_l2 ken\_l2, ken\_l3 ken\_l3, ken\_h1 ken\_h1, ken\_h2 ken\_h2) <br> |
|  subroutine | [**derlogden**](Derive__3d_8f90.md#function-derlogden) (logden logden, ld\_l1 ld\_l1, ld\_l2 ld\_l2, ld\_l3 ld\_l3, ld\_h1 ld\_h1, ld\_h2 ld\_h2, ld\_h3 ld\_h3) <br> |
|  subroutine | [**dermachnumber**](Derive__3d_8f90.md#function-dermachnumber) (mach mach, mach\_l1 mach\_l1, mach\_l2 mach\_l2, mach\_l3 mach\_l3, mach\_h1 mach\_h1) <br> |
|  subroutine | [**dermaggrav**](Derive__3d_8f90.md#function-dermaggrav) (maggrav maggrav, grav\_l1 grav\_l1, grav\_l2 grav\_l2, grav\_l3 grav\_l3, grav\_h1 grav\_h1) <br> |
|  subroutine | [**dermagmom**](Derive__3d_8f90.md#function-dermagmom) (magmom magmom, mom\_l1 mom\_l1, mom\_l2 mom\_l2, mom\_l3 mom\_l3, mom\_h1 mom\_h1, mom\_h2 mom\_h2) <br> |
|  subroutine | [**dermagvel**](Derive__3d_8f90.md#function-dermagvel) (magvel magvel, vel\_l1 vel\_l1, vel\_l2 vel\_l2, vel\_l3 vel\_l3, vel\_h1 vel\_h1, vel\_h2 vel\_h2) <br> |
|  subroutine | [**dermagvort**](Derive__3d_8f90.md#function-dermagvort) (vort vort, v\_l1 v\_l1, v\_l2 v\_l2, v\_l3 v\_l3, v\_h1 v\_h1, v\_h2 v\_h2, v\_h3 v\_h3, nv nv, dat dat, dat\_l1 dat\_l1, dat\_l2 dat\_l2, dat\_l3 dat\_l3, dat\_h1 dat\_h1, dat\_h2 dat\_h2) <br> |
|  subroutine | [**dermomt**](Derive__3d_8f90.md#function-dermomt) (vel vel, vel\_l1 vel\_l1, vel\_l2 vel\_l2, vel\_l3 vel\_l3, vel\_h1 vel\_h1, vel\_h2 vel\_h2, vel\_h3 vel\_h3) <br> |
|  subroutine | [**dernull**](Derive__3d_8f90.md#function-dernull) (kineng kineng, ken\_l1 ken\_l1, ken\_l2 ken\_l2, ken\_l3 ken\_l3, ken\_h1 ken\_h1, ken\_h2 ken\_h2) <br> |
|  subroutine | [**derpres**](Derive__3d_8f90.md#function-derpres) (p p, p\_l1 p\_l1, p\_l2 p\_l2, p\_l3 p\_l3, p\_h1 p\_h1, p\_h2 p\_h2, p\_h3 p\_h3, ncomp\_p ncomp\_p, u u, u\_l1 u\_l1, u\_l2 u\_l2, u\_l3 u\_l3, u\_h1 u\_h1, u\_h2 u\_h2, u\_h3 u\_h3, ncomp\_u ncomp\_u, lo lo, hi hi, domlo domlo, domhi domhi, dx dx, xlo xlo, time time, dt dt, bc bc, level level, grid\_no grid\_no) <br> |
|  subroutine | [**dersoundspeed**](Derive__3d_8f90.md#function-dersoundspeed) (c c, c\_l1 c\_l1, c\_l2 c\_l2, c\_l3 c\_l3, c\_h1 c\_h1, c\_h2 c\_h2, c\_h3 c\_h3, ncomp\_c ncomp\_c) <br> |
|  subroutine | [**derspec**](Derive__3d_8f90.md#function-derspec) (spec spec, spec\_l1 spec\_l1, spec\_l2 spec\_l2, spec\_l3 spec\_l3, spec\_h1 spec\_h1, spec\_h2 spec\_h2) <br> |
|  subroutine | [**derstate**](Derive__3d_8f90.md#function-derstate) (state state, state\_l1 state\_l1, state\_l2 state\_l2, state\_l3 state\_l3, state\_h1 state\_h1) <br> |
|  subroutine | [**dervel**](Derive__3d_8f90.md#function-dervel) (vel vel, vel\_l1 vel\_l1, vel\_l2 vel\_l2, vel\_l3 vel\_l3, vel\_h1 vel\_h1, vel\_h2 vel\_h2, vel\_h3 vel\_h3) <br> |








## Public Functions Documentation


### <a href="#function-ca-axangmom-x" id="function-ca-axangmom-x">function ca\_axangmom\_x </a>


```cpp
subroutine ca_axangmom_x (
    angmom_x angmom_x,
    angmom_x_l1 angmom_x_l1,
    angmom_x_l2 angmom_x_l2
) 
```



### <a href="#function-ca-axangmom-y" id="function-ca-axangmom-y">function ca\_axangmom\_y </a>


```cpp
subroutine ca_axangmom_y (
    angmom_y angmom_y,
    angmom_y_l1 angmom_y_l1,
    angmom_y_l2 angmom_y_l2
) 
```



### <a href="#function-ca-axangmom-z" id="function-ca-axangmom-z">function ca\_axangmom\_z </a>


```cpp
subroutine ca_axangmom_z (
    angmom_z angmom_z,
    angmom_z_l1 angmom_z_l1,
    angmom_z_l2 angmom_z_l2
) 
```



### <a href="#function-ca-axekin" id="function-ca-axekin">function ca\_axekin </a>


```cpp
subroutine ca_axekin (
    ekin ekin,
    ekin_l1 ekin_l1,
    ekin_l2 ekin_l2,
    ekin_l3 ekin_l3,
    ekin_h1 ekin_h1,
    ekin_h2 ekin_h2
) 
```



### <a href="#function-ca-axekinrho" id="function-ca-axekinrho">function ca\_axekinrho </a>


```cpp
subroutine ca_axekinrho (
    ekinrho ekinrho,
    ekinrho_l1 ekinrho_l1,
    ekinrho_l2 ekinrho_l2,
    ekinrho_l3 ekinrho_l3
) 
```



### <a href="#function-ca-axekinv" id="function-ca-axekinv">function ca\_axekinv </a>


```cpp
subroutine ca_axekinv (
    ekinv ekinv,
    ekinv_l1 ekinv_l1,
    ekinv_l2 ekinv_l2,
    ekinv_l3 ekinv_l3,
    ekinv_h1 ekinv_h1
) 
```



### <a href="#function-ca-axepot" id="function-ca-axepot">function ca\_axepot </a>


```cpp
subroutine ca_axepot (
    epot epot,
    epot_l1 epot_l1,
    epot_l2 epot_l2,
    epot_l3 epot_l3,
    epot_h1 epot_h1,
    epot_h2 epot_h2
) 
```



### <a href="#function-ca-axphase" id="function-ca-axphase">function ca\_axphase </a>


```cpp
subroutine ca_axphase (
    phase phase,
    phase_l1 phase_l1,
    phase_l2 phase_l2,
    phase_l3 phase_l3,
    phase_h1 phase_h1
) 
```



### <a href="#function-ca-axvel" id="function-ca-axvel">function ca\_axvel </a>


```cpp
subroutine ca_axvel (
    ekinv ekinv,
    ekinv_l1 ekinv_l1,
    ekinv_l2 ekinv_l2,
    ekinv_l3 ekinv_l3,
    ekinv_h1 ekinv_h1
) 
```



### <a href="#function-ca-dererrx" id="function-ca-dererrx">function ca\_dererrx </a>


```cpp
subroutine ca_dererrx (
    err err,
    err_l1 err_l1,
    err_l2 err_l2,
    err_l3 err_l3,
    err_h1 err_h1,
    err_h2 err_h2
) 
```



### <a href="#function-ca-dererry" id="function-ca-dererry">function ca\_dererry </a>


```cpp
subroutine ca_dererry (
    err err,
    err_l1 err_l1,
    err_l2 err_l2,
    err_l3 err_l3,
    err_h1 err_h1,
    err_h2 err_h2
) 
```



### <a href="#function-ca-dererrz" id="function-ca-dererrz">function ca\_dererrz </a>


```cpp
subroutine ca_dererrz (
    err err,
    err_l1 err_l1,
    err_l2 err_l2,
    err_l3 err_l3,
    err_h1 err_h1,
    err_h2 err_h2
) 
```



### <a href="#function-derdivu" id="function-derdivu">function derdivu </a>


```cpp
subroutine derdivu (
    divu divu,
    div_l1 div_l1,
    div_l2 div_l2,
    div_l3 div_l3,
    div_h1 div_h1,
    div_h2 div_h2,
    div_h3 div_h3
) 
```



### <a href="#function-dereint1" id="function-dereint1">function dereint1 </a>


```cpp
subroutine dereint1 (
    e e,
    e_l1 e_l1,
    e_l2 e_l2,
    e_l3 e_l3,
    e_h1 e_h1,
    e_h2 e_h2,
    e_h3 e_h3,
    ncomp_e ncomp_e,
    u u,
    u_l1 u_l1,
    u_l2 u_l2,
    u_l3 u_l3,
    u_h1 u_h1,
    u_h2 u_h2,
    u_h3 u_h3,
    ncomp_u ncomp_u,
    lo lo,
    hi hi,
    domlo domlo,
    domhi domhi,
    dx dx,
    xlo xlo,
    time time,
    dt dt,
    bc bc,
    level level,
    grid_no grid_no
) 
```



### <a href="#function-dereint2" id="function-dereint2">function dereint2 </a>


```cpp
subroutine dereint2 (
    e e,
    e_l1 e_l1,
    e_l2 e_l2,
    e_l3 e_l3,
    e_h1 e_h1,
    e_h2 e_h2,
    e_h3 e_h3,
    ncomp_e ncomp_e,
    u u,
    u_l1 u_l1,
    u_l2 u_l2,
    u_l3 u_l3,
    u_h1 u_h1,
    u_h2 u_h2,
    u_h3 u_h3,
    ncomp_u ncomp_u,
    lo lo,
    hi hi,
    domlo domlo,
    domhi domhi,
    dx dx,
    xlo xlo,
    time time,
    dt dt,
    bc bc,
    level level,
    grid_no grid_no
) 
```



### <a href="#function-derentropy" id="function-derentropy">function derentropy </a>


```cpp
subroutine derentropy (
    s s,
    s_l1 s_l1,
    s_l2 s_l2,
    s_l3 s_l3,
    s_h1 s_h1,
    s_h2 s_h2,
    s_h3 s_h3,
    ncomp_s ncomp_s,
    u u,
    u_l1 u_l1,
    u_l2 u_l2,
    u_l3 u_l3,
    u_h1 u_h1,
    u_h2 u_h2,
    u_h3 u_h3,
    ncomp_u ncomp_u
) 
```



### <a href="#function-derkineng" id="function-derkineng">function derkineng </a>


```cpp
subroutine derkineng (
    kineng kineng,
    ken_l1 ken_l1,
    ken_l2 ken_l2,
    ken_l3 ken_l3,
    ken_h1 ken_h1,
    ken_h2 ken_h2
) 
```



### <a href="#function-derlogden" id="function-derlogden">function derlogden </a>


```cpp
subroutine derlogden (
    logden logden,
    ld_l1 ld_l1,
    ld_l2 ld_l2,
    ld_l3 ld_l3,
    ld_h1 ld_h1,
    ld_h2 ld_h2,
    ld_h3 ld_h3
) 
```



### <a href="#function-dermachnumber" id="function-dermachnumber">function dermachnumber </a>


```cpp
subroutine dermachnumber (
    mach mach,
    mach_l1 mach_l1,
    mach_l2 mach_l2,
    mach_l3 mach_l3,
    mach_h1 mach_h1
) 
```



### <a href="#function-dermaggrav" id="function-dermaggrav">function dermaggrav </a>


```cpp
subroutine dermaggrav (
    maggrav maggrav,
    grav_l1 grav_l1,
    grav_l2 grav_l2,
    grav_l3 grav_l3,
    grav_h1 grav_h1
) 
```



### <a href="#function-dermagmom" id="function-dermagmom">function dermagmom </a>


```cpp
subroutine dermagmom (
    magmom magmom,
    mom_l1 mom_l1,
    mom_l2 mom_l2,
    mom_l3 mom_l3,
    mom_h1 mom_h1,
    mom_h2 mom_h2
) 
```



### <a href="#function-dermagvel" id="function-dermagvel">function dermagvel </a>


```cpp
subroutine dermagvel (
    magvel magvel,
    vel_l1 vel_l1,
    vel_l2 vel_l2,
    vel_l3 vel_l3,
    vel_h1 vel_h1,
    vel_h2 vel_h2
) 
```



### <a href="#function-dermagvort" id="function-dermagvort">function dermagvort </a>


```cpp
subroutine dermagvort (
    vort vort,
    v_l1 v_l1,
    v_l2 v_l2,
    v_l3 v_l3,
    v_h1 v_h1,
    v_h2 v_h2,
    v_h3 v_h3,
    nv nv,
    dat dat,
    dat_l1 dat_l1,
    dat_l2 dat_l2,
    dat_l3 dat_l3,
    dat_h1 dat_h1,
    dat_h2 dat_h2
) 
```



### <a href="#function-dermomt" id="function-dermomt">function dermomt </a>


```cpp
subroutine dermomt (
    vel vel,
    vel_l1 vel_l1,
    vel_l2 vel_l2,
    vel_l3 vel_l3,
    vel_h1 vel_h1,
    vel_h2 vel_h2,
    vel_h3 vel_h3
) 
```



### <a href="#function-dernull" id="function-dernull">function dernull </a>


```cpp
subroutine dernull (
    kineng kineng,
    ken_l1 ken_l1,
    ken_l2 ken_l2,
    ken_l3 ken_l3,
    ken_h1 ken_h1,
    ken_h2 ken_h2
) 
```



### <a href="#function-derpres" id="function-derpres">function derpres </a>


```cpp
subroutine derpres (
    p p,
    p_l1 p_l1,
    p_l2 p_l2,
    p_l3 p_l3,
    p_h1 p_h1,
    p_h2 p_h2,
    p_h3 p_h3,
    ncomp_p ncomp_p,
    u u,
    u_l1 u_l1,
    u_l2 u_l2,
    u_l3 u_l3,
    u_h1 u_h1,
    u_h2 u_h2,
    u_h3 u_h3,
    ncomp_u ncomp_u,
    lo lo,
    hi hi,
    domlo domlo,
    domhi domhi,
    dx dx,
    xlo xlo,
    time time,
    dt dt,
    bc bc,
    level level,
    grid_no grid_no
) 
```



### <a href="#function-dersoundspeed" id="function-dersoundspeed">function dersoundspeed </a>


```cpp
subroutine dersoundspeed (
    c c,
    c_l1 c_l1,
    c_l2 c_l2,
    c_l3 c_l3,
    c_h1 c_h1,
    c_h2 c_h2,
    c_h3 c_h3,
    ncomp_c ncomp_c
) 
```



### <a href="#function-derspec" id="function-derspec">function derspec </a>


```cpp
subroutine derspec (
    spec spec,
    spec_l1 spec_l1,
    spec_l2 spec_l2,
    spec_l3 spec_l3,
    spec_h1 spec_h1,
    spec_h2 spec_h2
) 
```



### <a href="#function-derstate" id="function-derstate">function derstate </a>


```cpp
subroutine derstate (
    state state,
    state_l1 state_l1,
    state_l2 state_l2,
    state_l3 state_l3,
    state_h1 state_h1
) 
```



### <a href="#function-dervel" id="function-dervel">function dervel </a>


```cpp
subroutine dervel (
    vel vel,
    vel_l1 vel_l1,
    vel_l2 vel_l2,
    vel_l3 vel_l3,
    vel_h1 vel_h1,
    vel_h2 vel_h2,
    vel_h3 vel_h3
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/DerivedQuantities/Derive_3d.f90`