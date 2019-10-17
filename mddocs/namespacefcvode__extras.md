
# Namespace fcvode\_extras


[**Class List**](annotated.md) **>** [**fcvode\_extras**](namespacefcvode__extras.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**fcvode\_wrapper**](namespacefcvode__extras.md#function-fcvode-wrapper) (dt dt, rho\_in rho\_in, T\_in T\_in, ne\_in ne\_in, e\_in e\_in, neq neq, cvmem cvmem, sunvec\_y sunvec\_y, yvec yvec, T\_out T\_out, ne\_out ne\_out, e\_out e\_out) <br> |
|  subroutine | [**fcvode\_wrapper\_vec**](namespacefcvode__extras.md#function-fcvode-wrapper-vec) (dt dt, rho\_in rho\_in, T\_in T\_in, ne\_in ne\_in, e\_in e\_in, neq neq, cvmem cvmem, sunvec\_y sunvec\_y, yvec yvec, T\_out T\_out, ne\_out ne\_out, e\_out e\_out) <br> |
|  integer(c\_int) function | [**rhsfn**](namespacefcvode__extras.md#function-rhsfn) (tn tn, sunvec\_y sunvec\_y, sunvec\_f sunvec\_f, user\_data user\_data) <br> |
|  integer(c\_int) function | [**rhsfn\_vec**](namespacefcvode__extras.md#function-rhsfn-vec) (tn tn, sunvec\_y sunvec\_y, sunvec\_f sunvec\_f, user\_data user\_data) <br> |








## Public Functions Documentation


### <a href="#function-fcvode-wrapper" id="function-fcvode-wrapper">function fcvode\_wrapper </a>


```cpp
subroutine fcvode_extras::fcvode_wrapper (
    dt dt,
    rho_in rho_in,
    T_in T_in,
    ne_in ne_in,
    e_in e_in,
    neq neq,
    cvmem cvmem,
    sunvec_y sunvec_y,
    yvec yvec,
    T_out T_out,
    ne_out ne_out,
    e_out e_out
) 
```



### <a href="#function-fcvode-wrapper-vec" id="function-fcvode-wrapper-vec">function fcvode\_wrapper\_vec </a>


```cpp
subroutine fcvode_extras::fcvode_wrapper_vec (
    dt dt,
    rho_in rho_in,
    T_in T_in,
    ne_in ne_in,
    e_in e_in,
    neq neq,
    cvmem cvmem,
    sunvec_y sunvec_y,
    yvec yvec,
    T_out T_out,
    ne_out ne_out,
    e_out e_out
) 
```



### <a href="#function-rhsfn" id="function-rhsfn">function rhsfn </a>


```cpp
integer(c_int) function fcvode_extras::rhsfn (
    tn tn,
    sunvec_y sunvec_y,
    sunvec_f sunvec_f,
    user_data user_data
) 
```



### <a href="#function-rhsfn-vec" id="function-rhsfn-vec">function rhsfn\_vec </a>


```cpp
integer(c_int) function fcvode_extras::rhsfn_vec (
    tn tn,
    sunvec_y sunvec_y,
    sunvec_f sunvec_f,
    user_data user_data
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/HeatCool/fcvode_extras.f90`