
# File integrate\_state\_vode\_3d.f90


[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**integrate\_state\_vode\_3d.f90**](integrate__state__vode__3d_8f90.md)

[Go to the source code of this file.](integrate__state__vode__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**integrate\_state\_vode**](integrate__state__vode__3d_8f90.md#function-integrate-state-vode) (lo lo, hi hi, state state, s\_l1 s\_l1, s\_l2 s\_l2, s\_l3 s\_l3, s\_h1 s\_h1, s\_h2 s\_h2, s\_h3 s\_h3, diag\_eos diag\_eos, d\_l1 d\_l1, d\_l2 d\_l2, d\_l3 d\_l3, d\_h1 d\_h1, d\_h2 d\_h2, d\_h3 d\_h3, a a, half\_dt half\_dt, min\_iter min\_iter, max\_iter max\_iter) <br> |
|  subroutine | [**vode\_wrapper**](integrate__state__vode__3d_8f90.md#function-vode-wrapper) (dt dt, rho\_in rho\_in, T\_in T\_in, ne\_in ne\_in, e\_in e\_in, T\_out T\_out, ne\_out ne\_out, e\_out e\_out) <br> |








## Public Functions Documentation


### <a href="#function-integrate-state-vode" id="function-integrate-state-vode">function integrate\_state\_vode </a>


```cpp
subroutine integrate_state_vode (
    lo lo,
    hi hi,
    state state,
    s_l1 s_l1,
    s_l2 s_l2,
    s_l3 s_l3,
    s_h1 s_h1,
    s_h2 s_h2,
    s_h3 s_h3,
    diag_eos diag_eos,
    d_l1 d_l1,
    d_l2 d_l2,
    d_l3 d_l3,
    d_h1 d_h1,
    d_h2 d_h2,
    d_h3 d_h3,
    a a,
    half_dt half_dt,
    min_iter min_iter,
    max_iter max_iter
) 
```



### <a href="#function-vode-wrapper" id="function-vode-wrapper">function vode\_wrapper </a>


```cpp
subroutine vode_wrapper (
    dt dt,
    rho_in rho_in,
    T_in T_in,
    ne_in ne_in,
    e_in e_in,
    T_out T_out,
    ne_out ne_out,
    e_out e_out
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/HeatCool/integrate_state_vode_3d.f90`