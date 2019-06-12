
# Namespace cvode\_extras


[**Class List**](annotated.md) **>** [**cvode\_extras**](namespacecvode__extras.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**ode\_eos\_finalize**](namespacecvode__extras.md#function-ode-eos-finalize) (e\_out e\_out, rpar rpar, num\_eq num\_eq) <br> |
|  subroutine | [**ode\_eos\_setup**](namespacecvode__extras.md#function-ode-eos-setup) (a a, half\_dt half\_dt) <br> |
|  integer(c\_int) function | [**rhsfnreal**](namespacecvode__extras.md#function-rhsfnreal) (tn tn, yvec yvec, fvec fvec, rpar rpar, neq neq) <br> |








## Public Functions Documentation


### <a href="#function-ode-eos-finalize" id="function-ode-eos-finalize">function ode\_eos\_finalize </a>


```cpp
subroutine cvode_extras::ode_eos_finalize (
    e_out e_out,
    rpar rpar,
    num_eq num_eq
) 
```



### <a href="#function-ode-eos-setup" id="function-ode-eos-setup">function ode\_eos\_setup </a>


```cpp
subroutine cvode_extras::ode_eos_setup (
    a a,
    half_dt half_dt
) 
```



### <a href="#function-rhsfnreal" id="function-rhsfnreal">function rhsfnreal </a>


```cpp
integer(c_int) function cvode_extras::rhsfnreal (
    tn tn,
    yvec yvec,
    fvec fvec,
    rpar rpar,
    neq neq
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/HeatCool/cvode_extras.f90`