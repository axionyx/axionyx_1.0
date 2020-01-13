
# Namespace ct\_upwind


[**Class List**](annotated.md) **>** [**ct\_upwind**](namespacect__upwind.md)















## Classes

| Type | Name |
| ---: | :--- |
| interface | [**checkisnan**](interfacect__upwind_1_1checkisnan.md) <br> |





## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**checkisnanmult**](namespacect__upwind.md#function-checkisnanmult) (uout uout, uout\_l1 uout\_l1, uout\_l2 uout\_l2, uout\_l3 uout\_l3, uout\_h1 uout\_h1, uout\_h2 uout\_h2, uout\_h3 uout\_h3, num num) <br> |
|  subroutine | [**checkisnans**](namespacect__upwind.md#function-checkisnans) (uout uout, uout\_l1 uout\_l1, uout\_l2 uout\_l2, uout\_l3 uout\_l3, uout\_h1 uout\_h1, uout\_h2 uout\_h2, uout\_h3 uout\_h3) <br> |
|  subroutine | [**checknegdens**](namespacect__upwind.md#function-checknegdens) (uout uout, uout\_l1 uout\_l1, uout\_l2 uout\_l2, uout\_l3 uout\_l3, uout\_h1 uout\_h1, uout\_h2 uout\_h2, uout\_h3 uout\_h3) <br> |
|  subroutine | [**constoprim**](namespacect__upwind.md#function-constoprim) (q q, u u, q\_l1 q\_l1, q\_l2 q\_l2, q\_l3 q\_l3, q\_h1 q\_h1, q\_h2 q\_h2, q\_h3 q\_h3) <br> |
|  subroutine | [**corner\_couple**](namespacect__upwind.md#function-corner-couple) (lo lo, hi hi, uL uL, uR uR, um um, up up, q\_l1 q\_l1, q\_l2 q\_l2, q\_l3 q\_l3, q\_h1 q\_h1, q\_h2 q\_h2, q\_h3 q\_h3, flxx flxx, flxx\_l1 flxx\_l1, flxx\_l2 flxx\_l2, flxx\_l3 flxx\_l3, flxx\_h1 flxx\_h1, flxx\_h2 flxx\_h2, flxx\_h3 flxx\_h3, flxy flxy, flxy\_l1 flxy\_l1, flxy\_l2 flxy\_l2, flxy\_l3 flxy\_l3, flxy\_h1 flxy\_h1, flxy\_h2 flxy\_h2, flxy\_h3 flxy\_h3, flxz flxz, flxz\_l1 flxz\_l1, flxz\_l2 flxz\_l2, flxz\_l3 flxz\_l3, flxz\_h1 flxz\_h1, flxz\_h2 flxz\_h2, flxz\_h3 flxz\_h3, dx dx, dy dy, dz dz, dt dt) <br> |
|  subroutine | [**corner\_couple\_mag**](namespacect__upwind.md#function-corner-couple-mag) (lo lo, hi hi, uL uL, uR uR, um um, up up, q\_l1 q\_l1, q\_l2 q\_l2, q\_l3 q\_l3, q\_h1 q\_h1, q\_h2 q\_h2, q\_h3 q\_h3, Ex Ex, ex\_l1 ex\_l1, ex\_l2 ex\_l2, ex\_l3 ex\_l3, ex\_h1 ex\_h1, ex\_h2 ex\_h2, ex\_h3 ex\_h3, Ey Ey, ey\_l1 ey\_l1, ey\_l2 ey\_l2, ey\_l3 ey\_l3, ey\_h1 ey\_h1, ey\_h2 ey\_h2, ey\_h3 ey\_h3, Ez Ez, ez\_l1 ez\_l1, ez\_l2 ez\_l2, ez\_l3 ez\_l3, ez\_h1 ez\_h1, ez\_h2 ez\_h2, ez\_h3 ez\_h3, dx dx, dy dy, dz dz, dt dt) <br> |
|  subroutine, public | [**corner\_transport**](namespacect__upwind.md#function-corner-transport) (q q, qm qm, qp qp, q\_l1 q\_l1, q\_l2 q\_l2, q\_l3 q\_l3, q\_h1 q\_h1, q\_h2 q\_h2, q\_h3 q\_h3, flxx flxx, flxx\_l1 flxx\_l1, flxx\_l2 flxx\_l2, flxx\_l3 flxx\_l3, flxx\_h1 flxx\_h1, flxx\_h2 flxx\_h2, flxx\_h3 flxx\_h3, flxy flxy, flxy\_l1 flxy\_l1, flxy\_l2 flxy\_l2, flxy\_l3 flxy\_l3, flxy\_h1 flxy\_h1, flxy\_h2 flxy\_h2, flxy\_h3 flxy\_h3, flxz flxz, flxz\_l1 flxz\_l1, flxz\_l2 flxz\_l2, flxz\_l3 flxz\_l3, flxz\_h1 flxz\_h1, flxz\_h2 flxz\_h2, flxz\_h3 flxz\_h3, Ex Ex, ex\_l1 ex\_l1, ex\_l2 ex\_l2, ex\_l3 ex\_l3, ex\_h1 ex\_h1, ex\_h2 ex\_h2, ex\_h3 ex\_h3, Ey Ey, ey\_l1 ey\_l1, ey\_l2 ey\_l2, ey\_l3 ey\_l3, ey\_h1 ey\_h1, ey\_h2 ey\_h2, ey\_h3 ey\_h3, Ez Ez, ez\_l1 ez\_l1, ez\_l2 ez\_l2, ez\_l3 ez\_l3, ez\_h1 ez\_h1, ez\_h2 ez\_h2, ez\_h3 ez\_h3, dx dx, dy dy, dz dz, dt dt) <br> |
|  subroutine | [**half\_step**](namespacect__upwind.md#function-half-step) (lo lo, hi hi, uL uL, uR uR, um um, up up, q\_l1 q\_l1, q\_l2 q\_l2, q\_l3 q\_l3, q\_h1 q\_h1, q\_h2 q\_h2, q\_h3 q\_h3, flxx flxx, flxx\_l1 flxx\_l1, flxx\_l2 flxx\_l2, flxx\_l3 flxx\_l3, flxx\_h1 flxx\_h1, flxx\_h2 flxx\_h2, flxx\_h3 flxx\_h3, flxy flxy, flxy\_l1 flxy\_l1, flxy\_l2 flxy\_l2, flxy\_l3 flxy\_l3, flxy\_h1 flxy\_h1, flxy\_h2 flxy\_h2, flxy\_h3 flxy\_h3, flxz flxz, flxz\_l1 flxz\_l1, flxz\_l2 flxz\_l2, flxz\_l3 flxz\_l3, flxz\_h1 flxz\_h1, flxz\_h2 flxz\_h2, flxz\_h3 flxz\_h3, dx dx, dy dy, dz dz, dt dt) <br> |
|  subroutine | [**half\_step\_mag**](namespacect__upwind.md#function-half-step-mag) (lo lo, hi hi, uL uL, uR uR, um um, up up, q\_l1 q\_l1, q\_l2 q\_l2, q\_l3 q\_l3, q\_h1 q\_h1, q\_h2 q\_h2, q\_h3 q\_h3, Ex Ex, ex\_l1 ex\_l1, ex\_l2 ex\_l2, ex\_l3 ex\_l3, ex\_h1 ex\_h1, ex\_h2 ex\_h2, ex\_h3 ex\_h3, Ey Ey, ey\_l1 ey\_l1, ey\_l2 ey\_l2, ey\_l3 ey\_l3, ey\_h1 ey\_h1, ey\_h2 ey\_h2, ey\_h3 ey\_h3, Ez Ez, ez\_l1 ez\_l1, ez\_l2 ez\_l2, ez\_l3 ez\_l3, ez\_h1 ez\_h1, ez\_h2 ez\_h2, ez\_h3 ez\_h3, dx dx, dy dy, dz dz, dt dt) <br> |
|  subroutine | [**prim\_half**](namespacect__upwind.md#function-prim-half) (q2D q2D, q q, q\_l1 q\_l1, q\_l2 q\_l2, q\_l3 q\_l3, q\_h1 q\_h1, q\_h2 q\_h2, q\_h3 q\_h3, flxx flxx, flxx\_l1 flxx\_l1, flxx\_l2 flxx\_l2, flxx\_l3 flxx\_l3, flxx\_h1 flxx\_h1, flxx\_h2 flxx\_h2, flxx\_h3 flxx\_h3, flxy flxy, flxy\_l1 flxy\_l1, flxy\_l2 flxy\_l2, flxy\_l3 flxy\_l3, flxy\_h1 flxy\_h1, flxy\_h2 flxy\_h2, flxy\_h3 flxy\_h3, flxz flxz, flxz\_l1 flxz\_l1, flxz\_l2 flxz\_l2, flxz\_l3 flxz\_l3, flxz\_h1 flxz\_h1, flxz\_h2 flxz\_h2, flxz\_h3 flxz\_h3, dx dx, dy dy, dz dz, dt dt) <br> |
|  subroutine | [**qflux**](namespacect__upwind.md#function-qflux) (qflx qflx, flx flx, q q) <br> |








## Public Functions Documentation


### <a href="#function-checkisnanmult" id="function-checkisnanmult">function checkisnanmult </a>


```cpp
subroutine ct_upwind::checkisnanmult (
    uout uout,
    uout_l1 uout_l1,
    uout_l2 uout_l2,
    uout_l3 uout_l3,
    uout_h1 uout_h1,
    uout_h2 uout_h2,
    uout_h3 uout_h3,
    num num
) 
```



### <a href="#function-checkisnans" id="function-checkisnans">function checkisnans </a>


```cpp
subroutine ct_upwind::checkisnans (
    uout uout,
    uout_l1 uout_l1,
    uout_l2 uout_l2,
    uout_l3 uout_l3,
    uout_h1 uout_h1,
    uout_h2 uout_h2,
    uout_h3 uout_h3
) 
```



### <a href="#function-checknegdens" id="function-checknegdens">function checknegdens </a>


```cpp
subroutine ct_upwind::checknegdens (
    uout uout,
    uout_l1 uout_l1,
    uout_l2 uout_l2,
    uout_l3 uout_l3,
    uout_h1 uout_h1,
    uout_h2 uout_h2,
    uout_h3 uout_h3
) 
```



### <a href="#function-constoprim" id="function-constoprim">function constoprim </a>


```cpp
subroutine ct_upwind::constoprim (
    q q,
    u u,
    q_l1 q_l1,
    q_l2 q_l2,
    q_l3 q_l3,
    q_h1 q_h1,
    q_h2 q_h2,
    q_h3 q_h3
) 
```



### <a href="#function-corner-couple" id="function-corner-couple">function corner\_couple </a>


```cpp
subroutine ct_upwind::corner_couple (
    lo lo,
    hi hi,
    uL uL,
    uR uR,
    um um,
    up up,
    q_l1 q_l1,
    q_l2 q_l2,
    q_l3 q_l3,
    q_h1 q_h1,
    q_h2 q_h2,
    q_h3 q_h3,
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
    dx dx,
    dy dy,
    dz dz,
    dt dt
) 
```



### <a href="#function-corner-couple-mag" id="function-corner-couple-mag">function corner\_couple\_mag </a>


```cpp
subroutine ct_upwind::corner_couple_mag (
    lo lo,
    hi hi,
    uL uL,
    uR uR,
    um um,
    up up,
    q_l1 q_l1,
    q_l2 q_l2,
    q_l3 q_l3,
    q_h1 q_h1,
    q_h2 q_h2,
    q_h3 q_h3,
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
    dx dx,
    dy dy,
    dz dz,
    dt dt
) 
```



### <a href="#function-corner-transport" id="function-corner-transport">function corner\_transport </a>


```cpp
subroutine, public ct_upwind::corner_transport (
    q q,
    qm qm,
    qp qp,
    q_l1 q_l1,
    q_l2 q_l2,
    q_l3 q_l3,
    q_h1 q_h1,
    q_h2 q_h2,
    q_h3 q_h3,
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
    dx dx,
    dy dy,
    dz dz,
    dt dt
) 
```



### <a href="#function-half-step" id="function-half-step">function half\_step </a>


```cpp
subroutine ct_upwind::half_step (
    lo lo,
    hi hi,
    uL uL,
    uR uR,
    um um,
    up up,
    q_l1 q_l1,
    q_l2 q_l2,
    q_l3 q_l3,
    q_h1 q_h1,
    q_h2 q_h2,
    q_h3 q_h3,
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
    dx dx,
    dy dy,
    dz dz,
    dt dt
) 
```



### <a href="#function-half-step-mag" id="function-half-step-mag">function half\_step\_mag </a>


```cpp
subroutine ct_upwind::half_step_mag (
    lo lo,
    hi hi,
    uL uL,
    uR uR,
    um um,
    up up,
    q_l1 q_l1,
    q_l2 q_l2,
    q_l3 q_l3,
    q_h1 q_h1,
    q_h2 q_h2,
    q_h3 q_h3,
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
    dx dx,
    dy dy,
    dz dz,
    dt dt
) 
```



### <a href="#function-prim-half" id="function-prim-half">function prim\_half </a>


```cpp
subroutine ct_upwind::prim_half (
    q2D q2D,
    q q,
    q_l1 q_l1,
    q_l2 q_l2,
    q_l3 q_l3,
    q_h1 q_h1,
    q_h2 q_h2,
    q_h3 q_h3,
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
    dx dx,
    dy dy,
    dz dz,
    dt dt
) 
```



### <a href="#function-qflux" id="function-qflux">function qflux </a>


```cpp
subroutine ct_upwind::qflux (
    qflx qflx,
    flx flx,
    q q
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/MHD/ct_upwind.f90`