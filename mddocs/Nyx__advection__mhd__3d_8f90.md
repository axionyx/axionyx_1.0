
# File Nyx\_advection\_mhd\_3d.f90


[**File List**](files.md) **>** [**MHD**](dir_a5db59ee0cc93a408ad0433ba32613c6.md) **>** [**Nyx\_advection\_mhd\_3d.f90**](Nyx__advection__mhd__3d_8f90.md)

[Go to the source code of this file.](Nyx__advection__mhd__3d_8f90_source.md)












## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**bl\_constants\_module**](namespacebl__constants__module.md) <br> |






## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**consup**](Nyx__advection__mhd__3d_8f90.md#function-consup) (uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1, uin\_h2 uin\_h2, uin\_h3 uin\_h3, uout uout, uout\_l1 uout\_l1, uout\_l2 uout\_l2, uout\_l3 uout\_l3, uout\_h1 uout\_h1, uout\_h2 uout\_h2, uout\_h3 uout\_h3, bx bx, bx\_l1 bx\_l1, bx\_l2 bx\_l2, bx\_l3 bx\_l3, bx\_h1 bx\_h1, bx\_h2 bx\_h2, bx\_h3 bx\_h3, by by, by\_l1 by\_l1, by\_l2 by\_l2, by\_l3 by\_l3, by\_h1 by\_h1, by\_h2 by\_h2, by\_h3 by\_h3, bz bz, bz\_l1 bz\_l1, bz\_l2 bz\_l2, bz\_l3 bz\_l3, bz\_h1 bz\_h1, bz\_h2 bz\_h2, bz\_h3 bz\_h3, hydro\_src hydro\_src, hsrc\_l1 hsrc\_l1, hsrc\_l2 hsrc\_l2, hsrc\_l3 hsrc\_l3, hsrc\_h1 hsrc\_h1, hsrc\_h2 hsrc\_h2, hsrc\_h3 hsrc\_h3, lo lo, hi hi, dt dt, a\_old a\_old, a\_new a\_new) <br> |
|  subroutine | [**ctoprim**](Nyx__advection__mhd__3d_8f90.md#function-ctoprim) (lo lo, hi hi, uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1, uin\_h2 uin\_h2, uin\_h3 uin\_h3, bx bx, bxin\_l1 bxin\_l1, bxin\_l2 bxin\_l2, bxin\_l3 bxin\_l3, bxin\_h1 bxin\_h1, bxin\_h2 bxin\_h2, bxin\_h3 bxin\_h3, by by, byin\_l1 byin\_l1, byin\_l2 byin\_l2, byin\_l3 byin\_l3, byin\_h1 byin\_h1, byin\_h2 byin\_h2, byin\_h3 byin\_h3, bz bz, bzin\_l1 bzin\_l1, bzin\_l2 bzin\_l2, bzin\_l3 bzin\_l3, bzin\_h1 bzin\_h1, bzin\_h2 bzin\_h2, bzin\_h3 bzin\_h3, q q, cx cx, cy cy, cz cz, csml csml, flatn flatn, q\_l1 q\_l1, q\_l2 q\_l2, q\_l3 q\_l3, q\_h1 q\_h1, q\_h2 q\_h2, q\_h3 q\_h3, src src, src\_l1 src\_l1, src\_l2 src\_l2, src\_l3 src\_l3, src\_h1 src\_h1, src\_h2 src\_h2, src\_h3 src\_h3, srcQ srcQ, srcq\_l1 srcq\_l1, srcq\_l2 srcq\_l2, srcq\_l3 srcq\_l3, srcq\_h1 srcq\_h1, srcq\_h2 srcq\_h2, srcq\_h3 srcq\_h3, grav grav, gv\_l1 gv\_l1, gv\_l2 gv\_l2, gv\_l3 gv\_l3, gv\_h1 gv\_h1, gv\_h2 gv\_h2, gv\_h3 gv\_h3, courno courno, dx dx, dy dy, dz dz, dt dt, ngp ngp, ngf ngf, a\_old a\_old, a\_new a\_new) <br> |
|  subroutine | [**flux\_combo**](Nyx__advection__mhd__3d_8f90.md#function-flux-combo) (uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1, uin\_h2 uin\_h2, uin\_h3 uin\_h3, hydro\_src hydro\_src, hsrc\_l1 hsrc\_l1, hsrc\_l2 hsrc\_l2, hsrc\_l3 hsrc\_l3, hsrc\_h1 hsrc\_h1, hsrc\_h2 hsrc\_h2, hsrc\_h3 hsrc\_h3, flux1 flux1, flux1\_l1 flux1\_l1, flux1\_l2 flux1\_l2, flux1\_l3 flux1\_l3, flux1\_h1 flux1\_h1, flux1\_h2 flux1\_h2, flux1\_h3 flux1\_h3, flux2 flux2, flux2\_l1 flux2\_l1, flux2\_l2 flux2\_l2, flux2\_l3 flux2\_l3, flux2\_h1 flux2\_h1, flux2\_h2 flux2\_h2, flux2\_h3 flux2\_h3, flux3 flux3, flux3\_l1 flux3\_l1, flux3\_l2 flux3\_l2, flux3\_l3 flux3\_l3, flux3\_h1 flux3\_h1, flux3\_h2 flux3\_h2, flux3\_h3 flux3\_h3, divu\_cc divu\_cc, d\_l1 d\_l1, d\_l2 d\_l2, d\_l3 d\_l3, d\_h1 d\_h1, d\_h2 d\_h2, d\_h3 d\_h3, lo lo, hi hi, dx dx, dy dy, dz dz, dt dt, a\_old a\_old, a\_new a\_new) <br> |
|  subroutine | [**flux\_ener\_corr**](Nyx__advection__mhd__3d_8f90.md#function-flux-ener-corr) (hydro\_src hydro\_src, hsrc\_l1 hsrc\_l1, hsrc\_l2 hsrc\_l2, hsrc\_l3 hsrc\_l3, hsrc\_h1 hsrc\_h1, hsrc\_h2 hsrc\_h2, hsrc\_h3 hsrc\_h3, bcc bcc, bcc\_l1 bcc\_l1, bcc\_l2 bcc\_l2, bcc\_l3 bcc\_l3, bcc\_h1 bcc\_h1, bcc\_h2 bcc\_h2, bcc\_h3 bcc\_h3, flxx flxx, flxx\_l1 flxx\_l1, flxx\_l2 flxx\_l2, flxx\_l3 flxx\_l3, flxx\_h1 flxx\_h1, flxx\_h2 flxx\_h2, flxx\_h3 flxx\_h3, flxy flxy, flxy\_l1 flxy\_l1, flxy\_l2 flxy\_l2, flxy\_l3 flxy\_l3, flxy\_h1 flxy\_h1, flxy\_h2 flxy\_h2, flxy\_h3 flxy\_h3, flxz flxz, flxz\_l1 flxz\_l1, flxz\_l2 flxz\_l2, flxz\_l3 flxz\_l3, flxz\_h1 flxz\_h1, flxz\_h2 flxz\_h2, flxz\_h3 flxz\_h3, bxout bxout, bxout\_l1 bxout\_l1, bxout\_l2 bxout\_l2, bxout\_l3 bxout\_l3, bxout\_h1 bxout\_h1, bxout\_h2 bxout\_h2, bxout\_h3 bxout\_h3, byout byout, byout\_l1 byout\_l1, byout\_l2 byout\_l2, byout\_l3 byout\_l3, byout\_h1 byout\_h1, byout\_h2 byout\_h2, byout\_h3 byout\_h3, bzout bzout, bzout\_l1 bzout\_l1, bzout\_l2 bzout\_l2, bzout\_l3 bzout\_l3, bzout\_h1 bzout\_h1, bzout\_h2 bzout\_h2, bzout\_h3 bzout\_h3, lo lo, hi hi, dx dx, dy dy, dz dz, dt dt, a\_old a\_old, a\_new a\_new) <br> |
|  subroutine | [**fort\_make\_mhd\_sources**](Nyx__advection__mhd__3d_8f90.md#function-fort-make-mhd-sources) (time time, lo lo, hi hi, uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1, uin\_h2 uin\_h2, uin\_h3 uin\_h3, bxin bxin, bxin\_l1 bxin\_l1, bxin\_l2 bxin\_l2, bxin\_l3 bxin\_l3, bxin\_h1 bxin\_h1, bxin\_h2 bxin\_h2, bxin\_h3 bxin\_h3, byin byin, byin\_l1 byin\_l1, byin\_l2 byin\_l2, byin\_l3 byin\_l3, byin\_h1 byin\_h1, byin\_h2 byin\_h2, byin\_h3 byin\_h3, bzin bzin, bzin\_l1 bzin\_l1, bzin\_l2 bzin\_l2, bzin\_l3 bzin\_l3, bzin\_h1 bzin\_h1, bzin\_h2 bzin\_h2, bzin\_h3 bzin\_h3, bxout bxout, bxout\_l1 bxout\_l1, bxout\_l2 bxout\_l2, bxout\_l3 bxout\_l3, bxout\_h1 bxout\_h1, bxout\_h2 bxout\_h2, bxout\_h3 bxout\_h3, byout byout, byout\_l1 byout\_l1, byout\_l2 byout\_l2, byout\_l3 byout\_l3, byout\_h1 byout\_h1, byout\_h2 byout\_h2, byout\_h3 byout\_h3, bzout bzout, bzout\_l1 bzout\_l1, bzout\_l2 bzout\_l2, bzout\_l3 bzout\_l3, bzout\_h1 bzout\_h1, bzout\_h2 bzout\_h2, bzout\_h3 bzout\_h3, ugdnvx ugdnvx, ugdnvx\_l1 ugdnvx\_l1, ugdnvx\_l2 ugdnvx\_l2, ugdnvx\_l3 ugdnvx\_l3, ugdnvx\_h1 ugdnvx\_h1, ugdnvx\_h2 ugdnvx\_h2, ugdnvx\_h3 ugdnvx\_h3, ugdnvy ugdnvy, ugdnvy\_l1 ugdnvy\_l1, ugdnvy\_l2 ugdnvy\_l2, ugdnvy\_l3 ugdnvy\_l3, ugdnvy\_h1 ugdnvy\_h1, ugdnvy\_h2 ugdnvy\_h2, ugdnvy\_h3 ugdnvy\_h3, ugdnvz ugdnvz, ugdnvz\_l1 ugdnvz\_l1, ugdnvz\_l2 ugdnvz\_l2, ugdnvz\_l3 ugdnvz\_l3, ugdnvz\_h1 ugdnvz\_h1, ugdnvz\_h2 ugdnvz\_h2, ugdnvz\_h3 ugdnvz\_h3, src src, src\_l1 src\_l1, src\_l2 src\_l2, src\_l3 src\_l3, src\_h1 src\_h1, src\_h2 src\_h2, src\_h3 src\_h3, hydro\_src hydro\_src, hsrc\_l1 hsrc\_l1, hsrc\_l2 hsrc\_l2, hsrc\_l3 hsrc\_l3, hsrc\_h1 hsrc\_h1, hsrc\_h2 hsrc\_h2, hsrc\_h3 hsrc\_h3, divu\_cc divu\_cc, d\_l1 d\_l1, d\_l2 d\_l2, d\_l3 d\_l3, d\_h1 d\_h1, d\_h2 d\_h2, d\_h3 d\_h3, grav grav, gv\_l1 gv\_l1, gv\_l2 gv\_l2, gv\_l3 gv\_l3, gv\_h1 gv\_h1, gv\_h2 gv\_h2, gv\_h3 gv\_h3, delta delta, dt dt, flux1 flux1, flux1\_l1 flux1\_l1, flux1\_l2 flux1\_l2, flux1\_l3 flux1\_l3, flux1\_h1 flux1\_h1, flux1\_h2 flux1\_h2, flux1\_h3 flux1\_h3, flux2 flux2, flux2\_l1 flux2\_l1, flux2\_l2 flux2\_l2, flux2\_l3 flux2\_l3, flux2\_h1 flux2\_h1, flux2\_h2 flux2\_h2, flux2\_h3 flux2\_h3, flux3 flux3, flux3\_l1 flux3\_l1, flux3\_l2 flux3\_l2, flux3\_l3 flux3\_l3, flux3\_h1 flux3\_h1, flux3\_h2 flux3\_h2, flux3\_h3 flux3\_h3, Ex Ex, ex\_l1 ex\_l1, ex\_l2 ex\_l2, ex\_l3 ex\_l3, ex\_h1 ex\_h1, ex\_h2 ex\_h2, ex\_h3 ex\_h3, Ey Ey, ey\_l1 ey\_l1, ey\_l2 ey\_l2, ey\_l3 ey\_l3, ey\_h1 ey\_h1, ey\_h2 ey\_h2, ey\_h3 ey\_h3, Ez Ez, ez\_l1 ez\_l1, ez\_l2 ez\_l2, ez\_l3 ez\_l3, ez\_h1 ez\_h1, ez\_h2 ez\_h2, ez\_h3 ez\_h3, courno courno, a\_old a\_old, a\_new a\_new, print\_fortran\_warnings print\_fortran\_warnings) <br> |
|  subroutine | [**fort\_update\_mhd\_state**](Nyx__advection__mhd__3d_8f90.md#function-fort-update-mhd-state) (lo lo, hi hi, uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1, uin\_h2 uin\_h2, uin\_h3 uin\_h3, uout uout, uout\_l1 uout\_l1, uout\_l2 uout\_l2, uout\_l3 uout\_l3, uout\_h1 uout\_h1, uout\_h2 uout\_h2, uout\_h3 uout\_h3, bxout bxout, bxout\_l1 bxout\_l1, bxout\_l2 bxout\_l2, bxout\_l3 bxout\_l3, bxout\_h1 bxout\_h1, bxout\_h2 bxout\_h2, bxout\_h3 bxout\_h3, byout byout, byout\_l1 byout\_l1, byout\_l2 byout\_l2, byout\_l3 byout\_l3, byout\_h1 byout\_h1, byout\_h2 byout\_h2, byout\_h3 byout\_h3, bzout bzout, bzout\_l1 bzout\_l1, bzout\_l2 bzout\_l2, bzout\_l3 bzout\_l3, bzout\_h1 bzout\_h1, bzout\_h2 bzout\_h2, bzout\_h3 bzout\_h3, src src, src\_l1 src\_l1, src\_l2 src\_l2, src\_l3 src\_l3, src\_h1 src\_h1, src\_h2 src\_h2, src\_h3 src\_h3, hydro\_src hydro\_src, hsrc\_l1 hsrc\_l1, hsrc\_l2 hsrc\_l2, hsrc\_l3 hsrc\_l3, hsrc\_h1 hsrc\_h1, hsrc\_h2 hsrc\_h2, hsrc\_h3 hsrc\_h3, divu\_cc divu\_cc, d\_l1 d\_l1, d\_l2 d\_l2, d\_l3 d\_l3, d\_h1 d\_h1, d\_h2 d\_h2, d\_h3 d\_h3, dt dt, dx dx, dy dy, dz dz, a\_old a\_old, a\_new a\_new, print\_fortran\_warnings print\_fortran\_warnings) <br> |
|  subroutine | [**gen\_bcc**](Nyx__advection__mhd__3d_8f90.md#function-gen-bcc) (lo lo, hi hi, bcc bcc, bcc\_l1 bcc\_l1, bcc\_l2 bcc\_l2, bcc\_l3 bcc\_l3, bcc\_h1 bcc\_h1, bcc\_h2 bcc\_h2, bcc\_h3 bcc\_h3, bxin bxin, bxin\_l1 bxin\_l1, bxin\_l2 bxin\_l2, bxin\_l3 bxin\_l3, bxin\_h1 bxin\_h1, bxin\_h2 bxin\_h2, bxin\_h3 bxin\_h3, byin byin, byin\_l1 byin\_l1, byin\_l2 byin\_l2, byin\_l3 byin\_l3, byin\_h1 byin\_h1, byin\_h2 byin\_h2, byin\_h3 byin\_h3, bzin bzin, bzin\_l1 bzin\_l1, bzin\_l2 bzin\_l2, bzin\_l3 bzin\_l3, bzin\_h1 bzin\_h1, bzin\_h2 bzin\_h2, bzin\_h3 bzin\_h3) <br> |
|  subroutine | [**magup**](Nyx__advection__mhd__3d_8f90.md#function-magup) (bxin bxin, bxin\_l1 bxin\_l1, bxin\_l2 bxin\_l2, bxin\_l3 bxin\_l3, bxin\_h1 bxin\_h1, bxin\_h2 bxin\_h2, bxin\_h3 bxin\_h3, byin byin, byin\_l1 byin\_l1, byin\_l2 byin\_l2, byin\_l3 byin\_l3, byin\_h1 byin\_h1, byin\_h2 byin\_h2, byin\_h3 byin\_h3, bzin bzin, bzin\_l1 bzin\_l1, bzin\_l2 bzin\_l2, bzin\_l3 bzin\_l3, bzin\_h1 bzin\_h1, bzin\_h2 bzin\_h2, bzin\_h3 bzin\_h3, bxout bxout, bxout\_l1 bxout\_l1, bxout\_l2 bxout\_l2, bxout\_l3 bxout\_l3, bxout\_h1 bxout\_h1, bxout\_h2 bxout\_h2, bxout\_h3 bxout\_h3, byout byout, byout\_l1 byout\_l1, byout\_l2 byout\_l2, byout\_l3 byout\_l3, byout\_h1 byout\_h1, byout\_h2 byout\_h2, byout\_h3 byout\_h3, bzout bzout, bzout\_l1 bzout\_l1, bzout\_l2 bzout\_l2, bzout\_l3 bzout\_l3, bzout\_h1 bzout\_h1, bzout\_h2 bzout\_h2, bzout\_h3 bzout\_h3, src src, src\_l1 src\_l1, src\_l2 src\_l2, src\_l3 src\_l3, src\_h1 src\_h1, src\_h2 src\_h2, src\_h3 src\_h3, Ex Ex, ex\_l1 ex\_l1, ex\_l2 ex\_l2, ex\_l3 ex\_l3, ex\_h1 ex\_h1, ex\_h2 ex\_h2, ex\_h3 ex\_h3, Ey Ey, ey\_l1 ey\_l1, ey\_l2 ey\_l2, ey\_l3 ey\_l3, ey\_h1 ey\_h1, ey\_h2 ey\_h2, ey\_h3 ey\_h3, Ez Ez, ez\_l1 ez\_l1, ez\_l2 ez\_l2, ez\_l3 ez\_l3, ez\_h1 ez\_h1, ez\_h2 ez\_h2, ez\_h3 ez\_h3, lo lo, hi hi, dx dx, dy dy, dz dz, dt dt, a\_old a\_old, a\_new a\_new) <br> |








## Public Functions Documentation


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
    uout uout,
    uout_l1 uout_l1,
    uout_l2 uout_l2,
    uout_l3 uout_l3,
    uout_h1 uout_h1,
    uout_h2 uout_h2,
    uout_h3 uout_h3,
    bx bx,
    bx_l1 bx_l1,
    bx_l2 bx_l2,
    bx_l3 bx_l3,
    bx_h1 bx_h1,
    bx_h2 bx_h2,
    bx_h3 bx_h3,
    by by,
    by_l1 by_l1,
    by_l2 by_l2,
    by_l3 by_l3,
    by_h1 by_h1,
    by_h2 by_h2,
    by_h3 by_h3,
    bz bz,
    bz_l1 bz_l1,
    bz_l2 bz_l2,
    bz_l3 bz_l3,
    bz_h1 bz_h1,
    bz_h2 bz_h2,
    bz_h3 bz_h3,
    hydro_src hydro_src,
    hsrc_l1 hsrc_l1,
    hsrc_l2 hsrc_l2,
    hsrc_l3 hsrc_l3,
    hsrc_h1 hsrc_h1,
    hsrc_h2 hsrc_h2,
    hsrc_h3 hsrc_h3,
    lo lo,
    hi hi,
    dt dt,
    a_old a_old,
    a_new a_new
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
    uin_l3 uin_l3,
    uin_h1 uin_h1,
    uin_h2 uin_h2,
    uin_h3 uin_h3,
    bx bx,
    bxin_l1 bxin_l1,
    bxin_l2 bxin_l2,
    bxin_l3 bxin_l3,
    bxin_h1 bxin_h1,
    bxin_h2 bxin_h2,
    bxin_h3 bxin_h3,
    by by,
    byin_l1 byin_l1,
    byin_l2 byin_l2,
    byin_l3 byin_l3,
    byin_h1 byin_h1,
    byin_h2 byin_h2,
    byin_h3 byin_h3,
    bz bz,
    bzin_l1 bzin_l1,
    bzin_l2 bzin_l2,
    bzin_l3 bzin_l3,
    bzin_h1 bzin_h1,
    bzin_h2 bzin_h2,
    bzin_h3 bzin_h3,
    q q,
    cx cx,
    cy cy,
    cz cz,
    csml csml,
    flatn flatn,
    q_l1 q_l1,
    q_l2 q_l2,
    q_l3 q_l3,
    q_h1 q_h1,
    q_h2 q_h2,
    q_h3 q_h3,
    src src,
    src_l1 src_l1,
    src_l2 src_l2,
    src_l3 src_l3,
    src_h1 src_h1,
    src_h2 src_h2,
    src_h3 src_h3,
    srcQ srcQ,
    srcq_l1 srcq_l1,
    srcq_l2 srcq_l2,
    srcq_l3 srcq_l3,
    srcq_h1 srcq_h1,
    srcq_h2 srcq_h2,
    srcq_h3 srcq_h3,
    grav grav,
    gv_l1 gv_l1,
    gv_l2 gv_l2,
    gv_l3 gv_l3,
    gv_h1 gv_h1,
    gv_h2 gv_h2,
    gv_h3 gv_h3,
    courno courno,
    dx dx,
    dy dy,
    dz dz,
    dt dt,
    ngp ngp,
    ngf ngf,
    a_old a_old,
    a_new a_new
) 
```



### <a href="#function-flux-combo" id="function-flux-combo">function flux\_combo </a>


```cpp
subroutine flux_combo (
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
    hsrc_h1 hsrc_h1,
    hsrc_h2 hsrc_h2,
    hsrc_h3 hsrc_h3,
    flux1 flux1,
    flux1_l1 flux1_l1,
    flux1_l2 flux1_l2,
    flux1_l3 flux1_l3,
    flux1_h1 flux1_h1,
    flux1_h2 flux1_h2,
    flux1_h3 flux1_h3,
    flux2 flux2,
    flux2_l1 flux2_l1,
    flux2_l2 flux2_l2,
    flux2_l3 flux2_l3,
    flux2_h1 flux2_h1,
    flux2_h2 flux2_h2,
    flux2_h3 flux2_h3,
    flux3 flux3,
    flux3_l1 flux3_l1,
    flux3_l2 flux3_l2,
    flux3_l3 flux3_l3,
    flux3_h1 flux3_h1,
    flux3_h2 flux3_h2,
    flux3_h3 flux3_h3,
    divu_cc divu_cc,
    d_l1 d_l1,
    d_l2 d_l2,
    d_l3 d_l3,
    d_h1 d_h1,
    d_h2 d_h2,
    d_h3 d_h3,
    lo lo,
    hi hi,
    dx dx,
    dy dy,
    dz dz,
    dt dt,
    a_old a_old,
    a_new a_new
) 
```



### <a href="#function-flux-ener-corr" id="function-flux-ener-corr">function flux\_ener\_corr </a>


```cpp
subroutine flux_ener_corr (
    hydro_src hydro_src,
    hsrc_l1 hsrc_l1,
    hsrc_l2 hsrc_l2,
    hsrc_l3 hsrc_l3,
    hsrc_h1 hsrc_h1,
    hsrc_h2 hsrc_h2,
    hsrc_h3 hsrc_h3,
    bcc bcc,
    bcc_l1 bcc_l1,
    bcc_l2 bcc_l2,
    bcc_l3 bcc_l3,
    bcc_h1 bcc_h1,
    bcc_h2 bcc_h2,
    bcc_h3 bcc_h3,
    flxx flxx,
    flxx_l1 flxx_l1,
    flxx_l2 flxx_l2,
    flxx_l3 flxx_l3,
    flxx_h1 flxx_h1,
    flxx_h2 flxx_h2,
    flxx_h3 flxx_h3,
    flxy flxy,
    flxy_l1 flxy_l1,
    flxy_l2 flxy_l2,
    flxy_l3 flxy_l3,
    flxy_h1 flxy_h1,
    flxy_h2 flxy_h2,
    flxy_h3 flxy_h3,
    flxz flxz,
    flxz_l1 flxz_l1,
    flxz_l2 flxz_l2,
    flxz_l3 flxz_l3,
    flxz_h1 flxz_h1,
    flxz_h2 flxz_h2,
    flxz_h3 flxz_h3,
    bxout bxout,
    bxout_l1 bxout_l1,
    bxout_l2 bxout_l2,
    bxout_l3 bxout_l3,
    bxout_h1 bxout_h1,
    bxout_h2 bxout_h2,
    bxout_h3 bxout_h3,
    byout byout,
    byout_l1 byout_l1,
    byout_l2 byout_l2,
    byout_l3 byout_l3,
    byout_h1 byout_h1,
    byout_h2 byout_h2,
    byout_h3 byout_h3,
    bzout bzout,
    bzout_l1 bzout_l1,
    bzout_l2 bzout_l2,
    bzout_l3 bzout_l3,
    bzout_h1 bzout_h1,
    bzout_h2 bzout_h2,
    bzout_h3 bzout_h3,
    lo lo,
    hi hi,
    dx dx,
    dy dy,
    dz dz,
    dt dt,
    a_old a_old,
    a_new a_new
) 
```



### <a href="#function-fort-make-mhd-sources" id="function-fort-make-mhd-sources">function fort\_make\_mhd\_sources </a>


```cpp
subroutine fort_make_mhd_sources (
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
    bxin bxin,
    bxin_l1 bxin_l1,
    bxin_l2 bxin_l2,
    bxin_l3 bxin_l3,
    bxin_h1 bxin_h1,
    bxin_h2 bxin_h2,
    bxin_h3 bxin_h3,
    byin byin,
    byin_l1 byin_l1,
    byin_l2 byin_l2,
    byin_l3 byin_l3,
    byin_h1 byin_h1,
    byin_h2 byin_h2,
    byin_h3 byin_h3,
    bzin bzin,
    bzin_l1 bzin_l1,
    bzin_l2 bzin_l2,
    bzin_l3 bzin_l3,
    bzin_h1 bzin_h1,
    bzin_h2 bzin_h2,
    bzin_h3 bzin_h3,
    bxout bxout,
    bxout_l1 bxout_l1,
    bxout_l2 bxout_l2,
    bxout_l3 bxout_l3,
    bxout_h1 bxout_h1,
    bxout_h2 bxout_h2,
    bxout_h3 bxout_h3,
    byout byout,
    byout_l1 byout_l1,
    byout_l2 byout_l2,
    byout_l3 byout_l3,
    byout_h1 byout_h1,
    byout_h2 byout_h2,
    byout_h3 byout_h3,
    bzout bzout,
    bzout_l1 bzout_l1,
    bzout_l2 bzout_l2,
    bzout_l3 bzout_l3,
    bzout_h1 bzout_h1,
    bzout_h2 bzout_h2,
    bzout_h3 bzout_h3,
    ugdnvx ugdnvx,
    ugdnvx_l1 ugdnvx_l1,
    ugdnvx_l2 ugdnvx_l2,
    ugdnvx_l3 ugdnvx_l3,
    ugdnvx_h1 ugdnvx_h1,
    ugdnvx_h2 ugdnvx_h2,
    ugdnvx_h3 ugdnvx_h3,
    ugdnvy ugdnvy,
    ugdnvy_l1 ugdnvy_l1,
    ugdnvy_l2 ugdnvy_l2,
    ugdnvy_l3 ugdnvy_l3,
    ugdnvy_h1 ugdnvy_h1,
    ugdnvy_h2 ugdnvy_h2,
    ugdnvy_h3 ugdnvy_h3,
    ugdnvz ugdnvz,
    ugdnvz_l1 ugdnvz_l1,
    ugdnvz_l2 ugdnvz_l2,
    ugdnvz_l3 ugdnvz_l3,
    ugdnvz_h1 ugdnvz_h1,
    ugdnvz_h2 ugdnvz_h2,
    ugdnvz_h3 ugdnvz_h3,
    src src,
    src_l1 src_l1,
    src_l2 src_l2,
    src_l3 src_l3,
    src_h1 src_h1,
    src_h2 src_h2,
    src_h3 src_h3,
    hydro_src hydro_src,
    hsrc_l1 hsrc_l1,
    hsrc_l2 hsrc_l2,
    hsrc_l3 hsrc_l3,
    hsrc_h1 hsrc_h1,
    hsrc_h2 hsrc_h2,
    hsrc_h3 hsrc_h3,
    divu_cc divu_cc,
    d_l1 d_l1,
    d_l2 d_l2,
    d_l3 d_l3,
    d_h1 d_h1,
    d_h2 d_h2,
    d_h3 d_h3,
    grav grav,
    gv_l1 gv_l1,
    gv_l2 gv_l2,
    gv_l3 gv_l3,
    gv_h1 gv_h1,
    gv_h2 gv_h2,
    gv_h3 gv_h3,
    delta delta,
    dt dt,
    flux1 flux1,
    flux1_l1 flux1_l1,
    flux1_l2 flux1_l2,
    flux1_l3 flux1_l3,
    flux1_h1 flux1_h1,
    flux1_h2 flux1_h2,
    flux1_h3 flux1_h3,
    flux2 flux2,
    flux2_l1 flux2_l1,
    flux2_l2 flux2_l2,
    flux2_l3 flux2_l3,
    flux2_h1 flux2_h1,
    flux2_h2 flux2_h2,
    flux2_h3 flux2_h3,
    flux3 flux3,
    flux3_l1 flux3_l1,
    flux3_l2 flux3_l2,
    flux3_l3 flux3_l3,
    flux3_h1 flux3_h1,
    flux3_h2 flux3_h2,
    flux3_h3 flux3_h3,
    Ex Ex,
    ex_l1 ex_l1,
    ex_l2 ex_l2,
    ex_l3 ex_l3,
    ex_h1 ex_h1,
    ex_h2 ex_h2,
    ex_h3 ex_h3,
    Ey Ey,
    ey_l1 ey_l1,
    ey_l2 ey_l2,
    ey_l3 ey_l3,
    ey_h1 ey_h1,
    ey_h2 ey_h2,
    ey_h3 ey_h3,
    Ez Ez,
    ez_l1 ez_l1,
    ez_l2 ez_l2,
    ez_l3 ez_l3,
    ez_h1 ez_h1,
    ez_h2 ez_h2,
    ez_h3 ez_h3,
    courno courno,
    a_old a_old,
    a_new a_new,
    print_fortran_warnings print_fortran_warnings
) 
```



### <a href="#function-fort-update-mhd-state" id="function-fort-update-mhd-state">function fort\_update\_mhd\_state </a>


```cpp
subroutine fort_update_mhd_state (
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
    bxout bxout,
    bxout_l1 bxout_l1,
    bxout_l2 bxout_l2,
    bxout_l3 bxout_l3,
    bxout_h1 bxout_h1,
    bxout_h2 bxout_h2,
    bxout_h3 bxout_h3,
    byout byout,
    byout_l1 byout_l1,
    byout_l2 byout_l2,
    byout_l3 byout_l3,
    byout_h1 byout_h1,
    byout_h2 byout_h2,
    byout_h3 byout_h3,
    bzout bzout,
    bzout_l1 bzout_l1,
    bzout_l2 bzout_l2,
    bzout_l3 bzout_l3,
    bzout_h1 bzout_h1,
    bzout_h2 bzout_h2,
    bzout_h3 bzout_h3,
    src src,
    src_l1 src_l1,
    src_l2 src_l2,
    src_l3 src_l3,
    src_h1 src_h1,
    src_h2 src_h2,
    src_h3 src_h3,
    hydro_src hydro_src,
    hsrc_l1 hsrc_l1,
    hsrc_l2 hsrc_l2,
    hsrc_l3 hsrc_l3,
    hsrc_h1 hsrc_h1,
    hsrc_h2 hsrc_h2,
    hsrc_h3 hsrc_h3,
    divu_cc divu_cc,
    d_l1 d_l1,
    d_l2 d_l2,
    d_l3 d_l3,
    d_h1 d_h1,
    d_h2 d_h2,
    d_h3 d_h3,
    dt dt,
    dx dx,
    dy dy,
    dz dz,
    a_old a_old,
    a_new a_new,
    print_fortran_warnings print_fortran_warnings
) 
```



### <a href="#function-gen-bcc" id="function-gen-bcc">function gen\_bcc </a>


```cpp
subroutine gen_bcc (
    lo lo,
    hi hi,
    bcc bcc,
    bcc_l1 bcc_l1,
    bcc_l2 bcc_l2,
    bcc_l3 bcc_l3,
    bcc_h1 bcc_h1,
    bcc_h2 bcc_h2,
    bcc_h3 bcc_h3,
    bxin bxin,
    bxin_l1 bxin_l1,
    bxin_l2 bxin_l2,
    bxin_l3 bxin_l3,
    bxin_h1 bxin_h1,
    bxin_h2 bxin_h2,
    bxin_h3 bxin_h3,
    byin byin,
    byin_l1 byin_l1,
    byin_l2 byin_l2,
    byin_l3 byin_l3,
    byin_h1 byin_h1,
    byin_h2 byin_h2,
    byin_h3 byin_h3,
    bzin bzin,
    bzin_l1 bzin_l1,
    bzin_l2 bzin_l2,
    bzin_l3 bzin_l3,
    bzin_h1 bzin_h1,
    bzin_h2 bzin_h2,
    bzin_h3 bzin_h3
) 
```



### <a href="#function-magup" id="function-magup">function magup </a>


```cpp
subroutine magup (
    bxin bxin,
    bxin_l1 bxin_l1,
    bxin_l2 bxin_l2,
    bxin_l3 bxin_l3,
    bxin_h1 bxin_h1,
    bxin_h2 bxin_h2,
    bxin_h3 bxin_h3,
    byin byin,
    byin_l1 byin_l1,
    byin_l2 byin_l2,
    byin_l3 byin_l3,
    byin_h1 byin_h1,
    byin_h2 byin_h2,
    byin_h3 byin_h3,
    bzin bzin,
    bzin_l1 bzin_l1,
    bzin_l2 bzin_l2,
    bzin_l3 bzin_l3,
    bzin_h1 bzin_h1,
    bzin_h2 bzin_h2,
    bzin_h3 bzin_h3,
    bxout bxout,
    bxout_l1 bxout_l1,
    bxout_l2 bxout_l2,
    bxout_l3 bxout_l3,
    bxout_h1 bxout_h1,
    bxout_h2 bxout_h2,
    bxout_h3 bxout_h3,
    byout byout,
    byout_l1 byout_l1,
    byout_l2 byout_l2,
    byout_l3 byout_l3,
    byout_h1 byout_h1,
    byout_h2 byout_h2,
    byout_h3 byout_h3,
    bzout bzout,
    bzout_l1 bzout_l1,
    bzout_l2 bzout_l2,
    bzout_l3 bzout_l3,
    bzout_h1 bzout_h1,
    bzout_h2 bzout_h2,
    bzout_h3 bzout_h3,
    src src,
    src_l1 src_l1,
    src_l2 src_l2,
    src_l3 src_l3,
    src_h1 src_h1,
    src_h2 src_h2,
    src_h3 src_h3,
    Ex Ex,
    ex_l1 ex_l1,
    ex_l2 ex_l2,
    ex_l3 ex_l3,
    ex_h1 ex_h1,
    ex_h2 ex_h2,
    ex_h3 ex_h3,
    Ey Ey,
    ey_l1 ey_l1,
    ey_l2 ey_l2,
    ey_l3 ey_l3,
    ey_h1 ey_h1,
    ey_h2 ey_h2,
    ey_h3 ey_h3,
    Ez Ez,
    ez_l1 ez_l1,
    ez_l2 ez_l2,
    ez_l3 ez_l3,
    ez_h1 ez_h1,
    ez_h2 ez_h2,
    ez_h3 ez_h3,
    lo lo,
    hi hi,
    dx dx,
    dy dy,
    dz dz,
    dt dt,
    a_old a_old,
    a_new a_new
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/MHD/Nyx_advection_mhd_3d.f90`