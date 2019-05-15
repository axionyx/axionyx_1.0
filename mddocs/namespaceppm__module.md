
# Namespace ppm\_module


[**Class List**](annotated.md) **>** [**ppm\_module**](namespaceppm__module.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**ppm**](namespaceppm__module.md#function-ppm) (s s, s\_l1 s\_l1, s\_l2 s\_l2, s\_l3 s\_l3, s\_h1 s\_h1, s\_h2 s\_h2, s\_h3 s\_h3, u u, cspd cspd, qd\_l1 qd\_l1, qd\_l2 qd\_l2, qd\_l3 qd\_l3, qd\_h1 qd\_h1, qd\_h2 qd\_h2, qd\_h3 qd\_h3, flatn flatn, f\_l1 f\_l1, f\_l2 f\_l2, f\_l3 f\_l3, f\_h1 f\_h1, f\_h2 f\_h2, f\_h3 f\_h3, Ip Ip, Im Im, ilo1 ilo1, ilo2 ilo2, ihi1 ihi1, ihi2 ihi2, dx dx, dy dy, dz dz, dt dt, k3d k3d, kc kc, a\_old a\_old) <br> |
|  subroutine | [**ppm\_type1**](namespaceppm__module.md#function-ppm-type1) (s s, s\_l1 s\_l1, s\_l2 s\_l2, s\_l3 s\_l3, s\_h1 s\_h1, s\_h2 s\_h2, s\_h3 s\_h3, u u, cspd cspd, qd\_l1 qd\_l1, qd\_l2 qd\_l2, qd\_l3 qd\_l3, qd\_h1 qd\_h1, qd\_h2 qd\_h2, qd\_h3 qd\_h3, flatn flatn, f\_l1 f\_l1, f\_l2 f\_l2, f\_l3 f\_l3, f\_h1 f\_h1, f\_h2 f\_h2, f\_h3 f\_h3, Ip Ip, Im Im, ilo1 ilo1, ilo2 ilo2, ihi1 ihi1, ihi2 ihi2, dx dx, dy dy, dz dz, dt\_over\_a dt\_over\_a, k3d k3d, kc kc) <br> |
|  subroutine | [**ppm\_type2**](namespaceppm__module.md#function-ppm-type2) (s s, s\_l1 s\_l1, s\_l2 s\_l2, s\_l3 s\_l3, s\_h1 s\_h1, s\_h2 s\_h2, s\_h3 s\_h3, u u, cspd cspd, qd\_l1 qd\_l1, qd\_l2 qd\_l2, qd\_l3 qd\_l3, qd\_h1 qd\_h1, qd\_h2 qd\_h2, qd\_h3 qd\_h3, flatn flatn, f\_l1 f\_l1, f\_l2 f\_l2, f\_l3 f\_l3, f\_h1 f\_h1, f\_h2 f\_h2, f\_h3 f\_h3, Ip Ip, Im Im, ilo1 ilo1, ilo2 ilo2, ihi1 ihi1, ihi2 ihi2, dx dx, dy dy, dz dz, dt\_over\_a dt\_over\_a, k3d k3d, kc kc) <br> |








## Public Functions Documentation


### <a href="#function-ppm" id="function-ppm">function ppm </a>


```cpp
subroutine, public ppm_module::ppm (
    s s,
    s_l1 s_l1,
    s_l2 s_l2,
    s_l3 s_l3,
    s_h1 s_h1,
    s_h2 s_h2,
    s_h3 s_h3,
    u u,
    cspd cspd,
    qd_l1 qd_l1,
    qd_l2 qd_l2,
    qd_l3 qd_l3,
    qd_h1 qd_h1,
    qd_h2 qd_h2,
    qd_h3 qd_h3,
    flatn flatn,
    f_l1 f_l1,
    f_l2 f_l2,
    f_l3 f_l3,
    f_h1 f_h1,
    f_h2 f_h2,
    f_h3 f_h3,
    Ip Ip,
    Im Im,
    ilo1 ilo1,
    ilo2 ilo2,
    ihi1 ihi1,
    ihi2 ihi2,
    dx dx,
    dy dy,
    dz dz,
    dt dt,
    k3d k3d,
    kc kc,
    a_old a_old
) 
```



### <a href="#function-ppm-type1" id="function-ppm-type1">function ppm\_type1 </a>


```cpp
subroutine ppm_module::ppm_type1 (
    s s,
    s_l1 s_l1,
    s_l2 s_l2,
    s_l3 s_l3,
    s_h1 s_h1,
    s_h2 s_h2,
    s_h3 s_h3,
    u u,
    cspd cspd,
    qd_l1 qd_l1,
    qd_l2 qd_l2,
    qd_l3 qd_l3,
    qd_h1 qd_h1,
    qd_h2 qd_h2,
    qd_h3 qd_h3,
    flatn flatn,
    f_l1 f_l1,
    f_l2 f_l2,
    f_l3 f_l3,
    f_h1 f_h1,
    f_h2 f_h2,
    f_h3 f_h3,
    Ip Ip,
    Im Im,
    ilo1 ilo1,
    ilo2 ilo2,
    ihi1 ihi1,
    ihi2 ihi2,
    dx dx,
    dy dy,
    dz dz,
    dt_over_a dt_over_a,
    k3d k3d,
    kc kc
) 
```



### <a href="#function-ppm-type2" id="function-ppm-type2">function ppm\_type2 </a>


```cpp
subroutine ppm_module::ppm_type2 (
    s s,
    s_l1 s_l1,
    s_l2 s_l2,
    s_l3 s_l3,
    s_h1 s_h1,
    s_h2 s_h2,
    s_h3 s_h3,
    u u,
    cspd cspd,
    qd_l1 qd_l1,
    qd_l2 qd_l2,
    qd_l3 qd_l3,
    qd_h1 qd_h1,
    qd_h2 qd_h2,
    qd_h3 qd_h3,
    flatn flatn,
    f_l1 f_l1,
    f_l2 f_l2,
    f_l3 f_l3,
    f_h1 f_h1,
    f_h2 f_h2,
    f_h3 f_h3,
    Ip Ip,
    Im Im,
    ilo1 ilo1,
    ilo2 ilo2,
    ihi1 ihi1,
    ihi2 ihi2,
    dx dx,
    dy dy,
    dz dz,
    dt_over_a dt_over_a,
    k3d k3d,
    kc kc
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/HydroFortran/ppm_3d.f90`