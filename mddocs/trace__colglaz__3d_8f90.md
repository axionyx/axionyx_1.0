
# File trace\_colglaz\_3d.f90


[**File List**](files.md) **>** [**HydroFortran**](dir_1fab266cd447ad3f3624320661f845f1.md) **>** [**trace\_colglaz\_3d.f90**](trace__colglaz__3d_8f90.md)

[Go to the source code of this file.](trace__colglaz__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**tracexy\_cg**](trace__colglaz__3d_8f90.md#function-tracexy-cg) (q q, c c, qd\_l1 qd\_l1, qd\_l2 qd\_l2, qd\_l3 qd\_l3, qd\_h1 qd\_h1, qd\_h2 qd\_h2, qd\_h3 qd\_h3, dqx dqx, dqy dqy, dq\_l1 dq\_l1, dq\_l2 dq\_l2, dq\_l3 dq\_l3, dq\_h1 dq\_h1, dq\_h2 dq\_h2, dq\_h3 dq\_h3, qxm qxm, qxp qxp, qym qym, qyp qyp, qpd\_l1 qpd\_l1, qpd\_l2 qpd\_l2, qpd\_l3 qpd\_l3, qpd\_h1 qpd\_h1, qpd\_h2 qpd\_h2, qpd\_h3 qpd\_h3, ilo1 ilo1, ilo2 ilo2, ihi1 ihi1, ihi2 ihi2, dx dx, dy dy, dt dt, kc kc, k3d k3d, a\_old a\_old) <br> |
|  subroutine | [**tracez\_cg**](trace__colglaz__3d_8f90.md#function-tracez-cg) (q q, c c, qd\_l1 qd\_l1, qd\_l2 qd\_l2, qd\_l3 qd\_l3, qd\_h1 qd\_h1, qd\_h2 qd\_h2, qd\_h3 qd\_h3, dqz dqz, dq\_l1 dq\_l1, dq\_l2 dq\_l2, dq\_l3 dq\_l3, dq\_h1 dq\_h1, dq\_h2 dq\_h2, dq\_h3 dq\_h3, qzm qzm, qzp qzp, qpd\_l1 qpd\_l1, qpd\_l2 qpd\_l2, qpd\_l3 qpd\_l3, qpd\_h1 qpd\_h1, qpd\_h2 qpd\_h2, qpd\_h3 qpd\_h3, ilo1 ilo1, ilo2 ilo2, ihi1 ihi1, ihi2 ihi2, dz dz, dt dt, km km, kc kc, k3d k3d, a\_old a\_old) <br> |








## Public Functions Documentation


### <a href="#function-tracexy-cg" id="function-tracexy-cg">function tracexy\_cg </a>


```cpp
subroutine tracexy_cg (
    q q,
    c c,
    qd_l1 qd_l1,
    qd_l2 qd_l2,
    qd_l3 qd_l3,
    qd_h1 qd_h1,
    qd_h2 qd_h2,
    qd_h3 qd_h3,
    dqx dqx,
    dqy dqy,
    dq_l1 dq_l1,
    dq_l2 dq_l2,
    dq_l3 dq_l3,
    dq_h1 dq_h1,
    dq_h2 dq_h2,
    dq_h3 dq_h3,
    qxm qxm,
    qxp qxp,
    qym qym,
    qyp qyp,
    qpd_l1 qpd_l1,
    qpd_l2 qpd_l2,
    qpd_l3 qpd_l3,
    qpd_h1 qpd_h1,
    qpd_h2 qpd_h2,
    qpd_h3 qpd_h3,
    ilo1 ilo1,
    ilo2 ilo2,
    ihi1 ihi1,
    ihi2 ihi2,
    dx dx,
    dy dy,
    dt dt,
    kc kc,
    k3d k3d,
    a_old a_old
) 
```



### <a href="#function-tracez-cg" id="function-tracez-cg">function tracez\_cg </a>


```cpp
subroutine tracez_cg (
    q q,
    c c,
    qd_l1 qd_l1,
    qd_l2 qd_l2,
    qd_l3 qd_l3,
    qd_h1 qd_h1,
    qd_h2 qd_h2,
    qd_h3 qd_h3,
    dqz dqz,
    dq_l1 dq_l1,
    dq_l2 dq_l2,
    dq_l3 dq_l3,
    dq_h1 dq_h1,
    dq_h2 dq_h2,
    dq_h3 dq_h3,
    qzm qzm,
    qzp qzp,
    qpd_l1 qpd_l1,
    qpd_l2 qpd_l2,
    qpd_l3 qpd_l3,
    qpd_h1 qpd_h1,
    qpd_h2 qpd_h2,
    qpd_h3 qpd_h3,
    ilo1 ilo1,
    ilo2 ilo2,
    ihi1 ihi1,
    ihi2 ihi2,
    dz dz,
    dt dt,
    km km,
    kc kc,
    k3d k3d,
    a_old a_old
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/HydroFortran/trace_colglaz_3d.f90`