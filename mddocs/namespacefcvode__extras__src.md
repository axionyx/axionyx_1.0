
# Namespace fcvode\_extras\_src


[**Class List**](annotated.md) **>** [**fcvode\_extras\_src**](namespacefcvode__extras__src.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**fcvode\_wrapper\_with\_source**](namespacefcvode__extras__src.md#function-fcvode-wrapper-with-source) (dt dt, rho\_in rho\_in, T\_in T\_in, ne\_in ne\_in, e\_in e\_in, neq neq, cvmem cvmem, sunvec\_y sunvec\_y, yvec yvec, rho\_out rho\_out, T\_out T\_out, ne\_out ne\_out, e\_out e\_out, rho\_src rho\_src, e\_src e\_src) <br> |
|  subroutine | [**fcvode\_wrapper\_with\_source\_single**](namespacefcvode__extras__src.md#function-fcvode-wrapper-with-source-single) (dt dt, rho\_in rho\_in, T\_in T\_in, ne\_in ne\_in, e\_in e\_in, neq neq, cvmem cvmem, sunvec\_y sunvec\_y, yvec yvec, rho\_out rho\_out, T\_out T\_out, ne\_out ne\_out, e\_out e\_out, rho\_src rho\_src, e\_src e\_src) <br> |
|  integer(c\_int) function | [**rhsfn\_src**](namespacefcvode__extras__src.md#function-rhsfn-src) (tn tn, sunvec\_y sunvec\_y, sunvec\_f sunvec\_f, user\_data user\_data) <br> |








## Public Functions Documentation


### <a href="#function-fcvode-wrapper-with-source" id="function-fcvode-wrapper-with-source">function fcvode\_wrapper\_with\_source </a>


```cpp
subroutine fcvode_extras_src::fcvode_wrapper_with_source (
    dt dt,
    rho_in rho_in,
    T_in T_in,
    ne_in ne_in,
    e_in e_in,
    neq neq,
    cvmem cvmem,
    sunvec_y sunvec_y,
    yvec yvec,
    rho_out rho_out,
    T_out T_out,
    ne_out ne_out,
    e_out e_out,
    rho_src rho_src,
    e_src e_src
) 
```



### <a href="#function-fcvode-wrapper-with-source-single" id="function-fcvode-wrapper-with-source-single">function fcvode\_wrapper\_with\_source\_single </a>


```cpp
subroutine fcvode_extras_src::fcvode_wrapper_with_source_single (
    dt dt,
    rho_in rho_in,
    T_in T_in,
    ne_in ne_in,
    e_in e_in,
    neq neq,
    cvmem cvmem,
    sunvec_y sunvec_y,
    yvec yvec,
    rho_out rho_out,
    T_out T_out,
    ne_out ne_out,
    e_out e_out,
    rho_src rho_src,
    e_src e_src
) 
```



### <a href="#function-rhsfn-src" id="function-rhsfn-src">function rhsfn\_src </a>


```cpp
integer(c_int) function fcvode_extras_src::rhsfn_src (
    tn tn,
    sunvec_y sunvec_y,
    sunvec_f sunvec_f,
    user_data user_data
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/HeatCool/fcvode_extras_src.f90`