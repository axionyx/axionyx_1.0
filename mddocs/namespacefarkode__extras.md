
# Namespace farkode\_extras


[**Class List**](annotated.md) **>** [**farkode\_extras**](namespacefarkode__extras.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**farkode\_wrapper**](namespacefarkode__extras.md#function-farkode-wrapper) (dt dt, rho\_in rho\_in, T\_in T\_in, ne\_in ne\_in, e\_in e\_in, neq neq, cvmem cvmem, sunvec\_y sunvec\_y, yvec yvec, T\_out T\_out, ne\_out ne\_out, e\_out e\_out) <br> |
|  integer(c\_int) function | [**rhsfnark**](namespacefarkode__extras.md#function-rhsfnark) (tn tn, sunvec\_y sunvec\_y, sunvec\_f sunvec\_f, user\_data user\_data) <br> |








## Public Functions Documentation


### <a href="#function-farkode-wrapper" id="function-farkode-wrapper">function farkode\_wrapper </a>


```cpp
subroutine farkode_extras::farkode_wrapper (
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



### <a href="#function-rhsfnark" id="function-rhsfnark">function rhsfnark </a>


```cpp
integer(c_int) function farkode_extras::rhsfnark (
    tn tn,
    sunvec_y sunvec_y,
    sunvec_f sunvec_f,
    user_data user_data
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/HeatCool/farkode_extras.f90`