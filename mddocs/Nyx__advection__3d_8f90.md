
# File Nyx\_advection\_3d.f90


[**File List**](files.md) **>** [**HydroFortran**](dir_1fab266cd447ad3f3624320661f845f1.md) **>** [**Nyx\_advection\_3d.f90**](Nyx__advection__3d_8f90.md)

[Go to the source code of this file.](Nyx__advection__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**cmpflx**](Nyx__advection__3d_8f90.md#function-cmpflx) (qm qm, qp qp, qpd\_l1 qpd\_l1, qpd\_l2 qpd\_l2, qpd\_l3 qpd\_l3, qpd\_h1 qpd\_h1, qpd\_h2 qpd\_h2, qpd\_h3 qpd\_h3) <br> |
|  subroutine | [**consup**](Nyx__advection__3d_8f90.md#function-consup) (uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1, uin\_h2 uin\_h2, uin\_h3 uin\_h3, hydro\_src hydro\_src, hsrc\_l1 hsrc\_l1, hsrc\_l2 hsrc\_l2, hsrc\_l3 hsrc\_l3, hsrc\_h1 hsrc\_h1) <br> |
|  subroutine | [**ctoprim**](Nyx__advection__3d_8f90.md#function-ctoprim) (lo lo, hi hi, uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3) <br> |
|  subroutine | [**make\_divu\_nd**](Nyx__advection__3d_8f90.md#function-make-divu-nd) (lo lo, hi hi, q q, q\_l1 q\_l1, q\_l2 q\_l2, q\_l3 q\_l3, q\_h1 q\_h1, q\_h2 q\_h2, q\_h3 q\_h3, dx dx) <br> |
|  subroutine | [**riemannus**](Nyx__advection__3d_8f90.md#function-riemannus) (ql ql, qr qr, qpd\_l1 qpd\_l1, qpd\_l2 qpd\_l2, qpd\_l3 qpd\_l3, qpd\_h1 qpd\_h1, qpd\_h2 qpd\_h2) <br> |
|  subroutine | [**umeth3d**](Nyx__advection__3d_8f90.md#function-umeth3d) (q q, c c, csml csml, flatn flatn, qd\_l1 qd\_l1, qd\_l2 qd\_l2, qd\_l3 qd\_l3, qd\_h1 qd\_h1) <br> |








## Public Functions Documentation


### <a href="#function-cmpflx" id="function-cmpflx">function cmpflx </a>


```cpp
subroutine cmpflx (
    qm qm,
    qp qp,
    qpd_l1 qpd_l1,
    qpd_l2 qpd_l2,
    qpd_l3 qpd_l3,
    qpd_h1 qpd_h1,
    qpd_h2 qpd_h2,
    qpd_h3 qpd_h3
) 
```



### <a href="#function-consup" id="function-consup">function consup </a>


```cpp
subroutine consup (
    uin uin,
    uin_l1 uin_l1,
    uin_l2 uin_l2,
    uin_l3 uin_l3,
    uin_h1 uin_h1,
    uin_h2 uin_h2,
    uin_h3 uin_h3,
    hydro_src hydro_src,
    hsrc_l1 hsrc_l1,
    hsrc_l2 hsrc_l2,
    hsrc_l3 hsrc_l3,
    hsrc_h1 hsrc_h1
) 
```



### <a href="#function-ctoprim" id="function-ctoprim">function ctoprim </a>


```cpp
subroutine ctoprim (
    lo lo,
    hi hi,
    uin uin,
    uin_l1 uin_l1,
    uin_l2 uin_l2,
    uin_l3 uin_l3
) 
```



### <a href="#function-make-divu-nd" id="function-make-divu-nd">function make\_divu\_nd </a>


```cpp
subroutine make_divu_nd (
    lo lo,
    hi hi,
    q q,
    q_l1 q_l1,
    q_l2 q_l2,
    q_l3 q_l3,
    q_h1 q_h1,
    q_h2 q_h2,
    q_h3 q_h3,
    dx dx
) 
```



### <a href="#function-riemannus" id="function-riemannus">function riemannus </a>


```cpp
subroutine riemannus (
    ql ql,
    qr qr,
    qpd_l1 qpd_l1,
    qpd_l2 qpd_l2,
    qpd_l3 qpd_l3,
    qpd_h1 qpd_h1,
    qpd_h2 qpd_h2
) 
```



### <a href="#function-umeth3d" id="function-umeth3d">function umeth3d </a>


```cpp
subroutine umeth3d (
    q q,
    c c,
    csml csml,
    flatn flatn,
    qd_l1 qd_l1,
    qd_l2 qd_l2,
    qd_l3 qd_l3,
    qd_h1 qd_h1
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/HydroFortran/Nyx_advection_3d.f90`