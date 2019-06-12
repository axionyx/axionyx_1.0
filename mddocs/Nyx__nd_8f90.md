
# File Nyx\_nd.f90


[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Nyx\_nd.f90**](Nyx__nd_8f90.md)

[Go to the source code of this file.](Nyx__nd_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**fort\_get\_aux\_names**](Nyx__nd_8f90.md#function-fort-get-aux-names) (aux\_names aux\_names, iaux iaux, len len) <br> |
|  subroutine | [**fort\_get\_method\_params**](Nyx__nd_8f90.md#function-fort-get-method-params) (nGrowHyp nGrowHyp) <br> |
|  subroutine | [**fort\_get\_num\_aux**](Nyx__nd_8f90.md#function-fort-get-num-aux) (naux\_out naux\_out) <br> |
|  subroutine | [**fort\_get\_num\_spec**](Nyx__nd_8f90.md#function-fort-get-num-spec) (nspec\_out nspec\_out) <br> |
|  subroutine | [**fort\_get\_spec\_names**](Nyx__nd_8f90.md#function-fort-get-spec-names) (spec\_names spec\_names, ispec ispec, len len) <br> |
|  subroutine | [**fort\_init\_zhi**](Nyx__nd_8f90.md#function-fort-init-zhi) (lo lo, hi hi, nd nd, diag\_eos diag\_eos, d\_l1 d\_l1, d\_l2 d\_l2, d\_l3 d\_l3, d\_h1 d\_h1, d\_h2 d\_h2, d\_h3 d\_h3, ratio ratio, zhi zhi, z\_l1 z\_l1, z\_l2 z\_l2, z\_l3 z\_l3, z\_h1 z\_h1, z\_h2 z\_h2, z\_h3 z\_h3) <br> |
|  subroutine | [**fort\_network\_init**](Nyx__nd_8f90.md#function-fort-network-init) () <br> |
|  subroutine | [**fort\_set\_eos\_params**](Nyx__nd_8f90.md#function-fort-set-eos-params) (h\_species\_in h\_species\_in, he\_species\_in he\_species\_in) <br> |
|  subroutine | [**fort\_set\_method\_params**](Nyx__nd_8f90.md#function-fort-set-method-params) (dm dm, numadv numadv, ndiag\_in ndiag\_in, do\_hydro do\_hydro, ppm\_type\_in ppm\_type\_in) <br> |
|  subroutine | [**fort\_set\_problem\_params**](Nyx__nd_8f90.md#function-fort-set-problem-params) (dm dm, physbc\_lo\_in physbc\_lo\_in, physbc\_hi\_in physbc\_hi\_in, Outflow\_in Outflow\_in, Symmetry\_in Symmetry\_in) <br> |
|  subroutine | [**fort\_set\_small\_values**](Nyx__nd_8f90.md#function-fort-set-small-values) (average\_dens average\_dens, average\_temp average\_temp, a a, small\_dens\_inout small\_dens\_inout) <br> |
|  subroutine | [**fort\_set\_xhydrogen**](Nyx__nd_8f90.md#function-fort-set-xhydrogen) (xhydrogen\_in xhydrogen\_in) <br> |
|  integer(c\_int) function | [**get\_comp\_e\_int**](Nyx__nd_8f90.md#function-get-comp-e-int) () <br> |
|  integer(c\_int) function | [**get\_comp\_temp**](Nyx__nd_8f90.md#function-get-comp-temp) () <br> |
|  integer(c\_int) function | [**get\_comp\_urho**](Nyx__nd_8f90.md#function-get-comp-urho) () <br> |








## Public Functions Documentation


### <a href="#function-fort-get-aux-names" id="function-fort-get-aux-names">function fort\_get\_aux\_names </a>


```cpp
subroutine fort_get_aux_names (
    aux_names aux_names,
    iaux iaux,
    len len
) 
```



### <a href="#function-fort-get-method-params" id="function-fort-get-method-params">function fort\_get\_method\_params </a>


```cpp
subroutine fort_get_method_params (
    nGrowHyp nGrowHyp
) 
```



### <a href="#function-fort-get-num-aux" id="function-fort-get-num-aux">function fort\_get\_num\_aux </a>


```cpp
subroutine fort_get_num_aux (
    naux_out naux_out
) 
```



### <a href="#function-fort-get-num-spec" id="function-fort-get-num-spec">function fort\_get\_num\_spec </a>


```cpp
subroutine fort_get_num_spec (
    nspec_out nspec_out
) 
```



### <a href="#function-fort-get-spec-names" id="function-fort-get-spec-names">function fort\_get\_spec\_names </a>


```cpp
subroutine fort_get_spec_names (
    spec_names spec_names,
    ispec ispec,
    len len
) 
```



### <a href="#function-fort-init-zhi" id="function-fort-init-zhi">function fort\_init\_zhi </a>


```cpp
subroutine fort_init_zhi (
    lo lo,
    hi hi,
    nd nd,
    diag_eos diag_eos,
    d_l1 d_l1,
    d_l2 d_l2,
    d_l3 d_l3,
    d_h1 d_h1,
    d_h2 d_h2,
    d_h3 d_h3,
    ratio ratio,
    zhi zhi,
    z_l1 z_l1,
    z_l2 z_l2,
    z_l3 z_l3,
    z_h1 z_h1,
    z_h2 z_h2,
    z_h3 z_h3
) 
```



### <a href="#function-fort-network-init" id="function-fort-network-init">function fort\_network\_init </a>


```cpp
subroutine fort_network_init () 
```



### <a href="#function-fort-set-eos-params" id="function-fort-set-eos-params">function fort\_set\_eos\_params </a>


```cpp
subroutine fort_set_eos_params (
    h_species_in h_species_in,
    he_species_in he_species_in
) 
```



### <a href="#function-fort-set-method-params" id="function-fort-set-method-params">function fort\_set\_method\_params </a>


```cpp
subroutine fort_set_method_params (
    dm dm,
    numadv numadv,
    ndiag_in ndiag_in,
    do_hydro do_hydro,
    ppm_type_in ppm_type_in
) 
```



### <a href="#function-fort-set-problem-params" id="function-fort-set-problem-params">function fort\_set\_problem\_params </a>


```cpp
subroutine fort_set_problem_params (
    dm dm,
    physbc_lo_in physbc_lo_in,
    physbc_hi_in physbc_hi_in,
    Outflow_in Outflow_in,
    Symmetry_in Symmetry_in
) 
```



### <a href="#function-fort-set-small-values" id="function-fort-set-small-values">function fort\_set\_small\_values </a>


```cpp
subroutine fort_set_small_values (
    average_dens average_dens,
    average_temp average_temp,
    a a,
    small_dens_inout small_dens_inout
) 
```



### <a href="#function-fort-set-xhydrogen" id="function-fort-set-xhydrogen">function fort\_set\_xhydrogen </a>


```cpp
subroutine fort_set_xhydrogen (
    xhydrogen_in xhydrogen_in
) 
```



### <a href="#function-get-comp-e-int" id="function-get-comp-e-int">function get\_comp\_e\_int </a>


```cpp
integer(c_int) function get_comp_e_int () 
```



### <a href="#function-get-comp-temp" id="function-get-comp-temp">function get\_comp\_temp </a>


```cpp
integer(c_int) function get_comp_temp () 
```



### <a href="#function-get-comp-urho" id="function-get-comp-urho">function get\_comp\_urho </a>


```cpp
integer(c_int) function get_comp_urho () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/Nyx_nd.f90`