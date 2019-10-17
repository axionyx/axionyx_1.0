
# Namespace trace\_src\_module


[**Class List**](annotated.md) **>** [**trace\_src\_module**](namespacetrace__src__module.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**tracex\_src**](namespacetrace__src__module.md#function-tracex-src) (q q, c c, qd\_l1 qd\_l1, qd\_l2 qd\_l2, qd\_l3 qd\_l3, qd\_h1 qd\_h1, qd\_h2 qd\_h2, qd\_h3 qd\_h3, qxm qxm, qxp qxp, qpd\_l1 qpd\_l1, qpd\_l2 qpd\_l2, qpd\_l3 qpd\_l3, qpd\_h1 qpd\_h1, qpd\_h2 qpd\_h2, qpd\_h3 qpd\_h3, srcQ srcQ, srcq\_l1 srcq\_l1, srcq\_l2 srcq\_l2, srcq\_l3 srcq\_l3, srcq\_h1 srcq\_h1, srcq\_h2 srcq\_h2, srcq\_h3 srcq\_h3, ilo1 ilo1, ilo2 ilo2, ihi1 ihi1, ihi2 ihi2, dt dt, a\_old a\_old, kc kc, k3d k3d) <br> |
|  subroutine, public | [**tracey\_src**](namespacetrace__src__module.md#function-tracey-src) (q q, c c, qd\_l1 qd\_l1, qd\_l2 qd\_l2, qd\_l3 qd\_l3, qd\_h1 qd\_h1, qd\_h2 qd\_h2, qd\_h3 qd\_h3, qym qym, qyp qyp, qpd\_l1 qpd\_l1, qpd\_l2 qpd\_l2, qpd\_l3 qpd\_l3, qpd\_h1 qpd\_h1, qpd\_h2 qpd\_h2, qpd\_h3 qpd\_h3, srcQ srcQ, srcq\_l1 srcq\_l1, srcq\_l2 srcq\_l2, srcq\_l3 srcq\_l3, srcq\_h1 srcq\_h1, srcq\_h2 srcq\_h2, srcq\_h3 srcq\_h3, ilo1 ilo1, ilo2 ilo2, ihi1 ihi1, ihi2 ihi2, dt dt, a\_old a\_old, kc kc, k3d k3d) <br> |
|  subroutine, public | [**tracez\_src**](namespacetrace__src__module.md#function-tracez-src) (q q, c c, qd\_l1 qd\_l1, qd\_l2 qd\_l2, qd\_l3 qd\_l3, qd\_h1 qd\_h1, qd\_h2 qd\_h2, qd\_h3 qd\_h3, qzm qzm, qzp qzp, qpd\_l1 qpd\_l1, qpd\_l2 qpd\_l2, qpd\_l3 qpd\_l3, qpd\_h1 qpd\_h1, qpd\_h2 qpd\_h2, qpd\_h3 qpd\_h3, srcQ srcQ, srcq\_l1 srcq\_l1, srcq\_l2 srcq\_l2, srcq\_l3 srcq\_l3, srcq\_h1 srcq\_h1, srcq\_h2 srcq\_h2, srcq\_h3 srcq\_h3, ilo1 ilo1, ilo2 ilo2, ihi1 ihi1, ihi2 ihi2, dt dt, a\_old a\_old, km km, kc kc, k3d k3d) <br> |








## Public Functions Documentation


### <a href="#function-tracex-src" id="function-tracex-src">function tracex\_src </a>


```cpp
subroutine, public trace_src_module::tracex_src (
    q q,
    c c,
    qd_l1 qd_l1,
    qd_l2 qd_l2,
    qd_l3 qd_l3,
    qd_h1 qd_h1,
    qd_h2 qd_h2,
    qd_h3 qd_h3,
    qxm qxm,
    qxp qxp,
    qpd_l1 qpd_l1,
    qpd_l2 qpd_l2,
    qpd_l3 qpd_l3,
    qpd_h1 qpd_h1,
    qpd_h2 qpd_h2,
    qpd_h3 qpd_h3,
    srcQ srcQ,
    srcq_l1 srcq_l1,
    srcq_l2 srcq_l2,
    srcq_l3 srcq_l3,
    srcq_h1 srcq_h1,
    srcq_h2 srcq_h2,
    srcq_h3 srcq_h3,
    ilo1 ilo1,
    ilo2 ilo2,
    ihi1 ihi1,
    ihi2 ihi2,
    dt dt,
    a_old a_old,
    kc kc,
    k3d k3d
) 
```



### <a href="#function-tracey-src" id="function-tracey-src">function tracey\_src </a>


```cpp
subroutine, public trace_src_module::tracey_src (
    q q,
    c c,
    qd_l1 qd_l1,
    qd_l2 qd_l2,
    qd_l3 qd_l3,
    qd_h1 qd_h1,
    qd_h2 qd_h2,
    qd_h3 qd_h3,
    qym qym,
    qyp qyp,
    qpd_l1 qpd_l1,
    qpd_l2 qpd_l2,
    qpd_l3 qpd_l3,
    qpd_h1 qpd_h1,
    qpd_h2 qpd_h2,
    qpd_h3 qpd_h3,
    srcQ srcQ,
    srcq_l1 srcq_l1,
    srcq_l2 srcq_l2,
    srcq_l3 srcq_l3,
    srcq_h1 srcq_h1,
    srcq_h2 srcq_h2,
    srcq_h3 srcq_h3,
    ilo1 ilo1,
    ilo2 ilo2,
    ihi1 ihi1,
    ihi2 ihi2,
    dt dt,
    a_old a_old,
    kc kc,
    k3d k3d
) 
```



### <a href="#function-tracez-src" id="function-tracez-src">function tracez\_src </a>


```cpp
subroutine, public trace_src_module::tracez_src (
    q q,
    c c,
    qd_l1 qd_l1,
    qd_l2 qd_l2,
    qd_l3 qd_l3,
    qd_h1 qd_h1,
    qd_h2 qd_h2,
    qd_h3 qd_h3,
    qzm qzm,
    qzp qzp,
    qpd_l1 qpd_l1,
    qpd_l2 qpd_l2,
    qpd_l3 qpd_l3,
    qpd_h1 qpd_h1,
    qpd_h2 qpd_h2,
    qpd_h3 qpd_h3,
    srcQ srcQ,
    srcq_l1 srcq_l1,
    srcq_l2 srcq_l2,
    srcq_l3 srcq_l3,
    srcq_h1 srcq_h1,
    srcq_h2 srcq_h2,
    srcq_h3 srcq_h3,
    ilo1 ilo1,
    ilo2 ilo2,
    ihi1 ihi1,
    ihi2 ihi2,
    dt dt,
    a_old a_old,
    km km,
    kc kc,
    k3d k3d
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/HydroFortran/trace_src_3d.f90`