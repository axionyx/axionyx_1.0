
# Namespace comoving\_nd\_module


[**Class List**](annotated.md) **>** [**comoving\_nd\_module**](namespacecomoving__nd__module.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**enforce\_final\_a**](namespacecomoving__nd__module.md#function-enforce-final-a) (old\_a old\_a, new\_a new\_a, dt dt, final\_a final\_a) <br> |
|  subroutine | [**enforce\_percent\_change**](namespacecomoving__nd__module.md#function-enforce-percent-change) (old\_a old\_a, new\_a new\_a, dt dt, change\_allowed change\_allowed) <br> |
|  subroutine | [**fort\_est\_lindt\_comoving\_a**](namespacecomoving__nd__module.md#function-fort-est-lindt-comoving-a) (old\_a old\_a, new\_a new\_a, dt dt) <br> |
|  subroutine | [**fort\_est\_maxdt\_comoving\_a**](namespacecomoving__nd__module.md#function-fort-est-maxdt-comoving-a) (old\_a old\_a, dt dt) <br> |
|  subroutine | [**fort\_estdt\_comoving\_a**](namespacecomoving__nd__module.md#function-fort-estdt-comoving-a) (old\_a old\_a, new\_a new\_a, dt dt, change\_allowed change\_allowed, fixed\_da fixed\_da, final\_a final\_a, dt\_modified dt\_modified) <br> |
|  subroutine | [**fort\_get\_hubble**](namespacecomoving__nd__module.md#function-fort-get-hubble) (hubble hubble) <br> |
|  subroutine | [**fort\_integrate\_comoving\_a**](namespacecomoving__nd__module.md#function-fort-integrate-comoving-a) (old\_a old\_a, new\_a new\_a, dt dt) <br> |
|  subroutine | [**fort\_integrate\_comoving\_a\_to\_a**](namespacecomoving__nd__module.md#function-fort-integrate-comoving-a-to-a) (old\_a old\_a, a\_value a\_value, dt dt) <br> |
|  subroutine | [**fort\_integrate\_comoving\_a\_to\_z**](namespacecomoving__nd__module.md#function-fort-integrate-comoving-a-to-z) (old\_a old\_a, z\_value z\_value, dt dt) <br> |
|  subroutine | [**fort\_integrate\_time\_given\_a**](namespacecomoving__nd__module.md#function-fort-integrate-time-given-a) (a0 a0, a1 a1, dt dt) <br> |
|  subroutine | [**fort\_set\_hubble**](namespacecomoving__nd__module.md#function-fort-set-hubble) (hubble hubble) <br> |
|  subroutine | [**fort\_set\_omb**](namespacecomoving__nd__module.md#function-fort-set-omb) (omb omb) <br> |
|  subroutine | [**fort\_set\_omm**](namespacecomoving__nd__module.md#function-fort-set-omm) (omm omm) <br> |
|  real(rt) function | [**invez**](namespacecomoving__nd__module.md#function-invez) (H0 H0, Om Om, a a) <br> |








## Public Functions Documentation


### <a href="#function-enforce-final-a" id="function-enforce-final-a">function enforce\_final\_a </a>


```cpp
subroutine comoving_nd_module::enforce_final_a (
    old_a old_a,
    new_a new_a,
    dt dt,
    final_a final_a
) 
```



### <a href="#function-enforce-percent-change" id="function-enforce-percent-change">function enforce\_percent\_change </a>


```cpp
subroutine comoving_nd_module::enforce_percent_change (
    old_a old_a,
    new_a new_a,
    dt dt,
    change_allowed change_allowed
) 
```



### <a href="#function-fort-est-lindt-comoving-a" id="function-fort-est-lindt-comoving-a">function fort\_est\_lindt\_comoving\_a </a>


```cpp
subroutine comoving_nd_module::fort_est_lindt_comoving_a (
    old_a old_a,
    new_a new_a,
    dt dt
) 
```



### <a href="#function-fort-est-maxdt-comoving-a" id="function-fort-est-maxdt-comoving-a">function fort\_est\_maxdt\_comoving\_a </a>


```cpp
subroutine comoving_nd_module::fort_est_maxdt_comoving_a (
    old_a old_a,
    dt dt
) 
```



### <a href="#function-fort-estdt-comoving-a" id="function-fort-estdt-comoving-a">function fort\_estdt\_comoving\_a </a>


```cpp
subroutine comoving_nd_module::fort_estdt_comoving_a (
    old_a old_a,
    new_a new_a,
    dt dt,
    change_allowed change_allowed,
    fixed_da fixed_da,
    final_a final_a,
    dt_modified dt_modified
) 
```



### <a href="#function-fort-get-hubble" id="function-fort-get-hubble">function fort\_get\_hubble </a>


```cpp
subroutine comoving_nd_module::fort_get_hubble (
    hubble hubble
) 
```



### <a href="#function-fort-integrate-comoving-a" id="function-fort-integrate-comoving-a">function fort\_integrate\_comoving\_a </a>


```cpp
subroutine comoving_nd_module::fort_integrate_comoving_a (
    old_a old_a,
    new_a new_a,
    dt dt
) 
```



### <a href="#function-fort-integrate-comoving-a-to-a" id="function-fort-integrate-comoving-a-to-a">function fort\_integrate\_comoving\_a\_to\_a </a>


```cpp
subroutine comoving_nd_module::fort_integrate_comoving_a_to_a (
    old_a old_a,
    a_value a_value,
    dt dt
) 
```



### <a href="#function-fort-integrate-comoving-a-to-z" id="function-fort-integrate-comoving-a-to-z">function fort\_integrate\_comoving\_a\_to\_z </a>


```cpp
subroutine comoving_nd_module::fort_integrate_comoving_a_to_z (
    old_a old_a,
    z_value z_value,
    dt dt
) 
```



### <a href="#function-fort-integrate-time-given-a" id="function-fort-integrate-time-given-a">function fort\_integrate\_time\_given\_a </a>


```cpp
subroutine comoving_nd_module::fort_integrate_time_given_a (
    a0 a0,
    a1 a1,
    dt dt
) 
```



### <a href="#function-fort-set-hubble" id="function-fort-set-hubble">function fort\_set\_hubble </a>


```cpp
subroutine comoving_nd_module::fort_set_hubble (
    hubble hubble
) 
```



### <a href="#function-fort-set-omb" id="function-fort-set-omb">function fort\_set\_omb </a>


```cpp
subroutine comoving_nd_module::fort_set_omb (
    omb omb
) 
```



### <a href="#function-fort-set-omm" id="function-fort-set-omm">function fort\_set\_omm </a>


```cpp
subroutine comoving_nd_module::fort_set_omm (
    omm omm
) 
```



### <a href="#function-invez" id="function-invez">function invez </a>


```cpp
real(rt) function comoving_nd_module::invez (
    H0 H0,
    Om Om,
    a a
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/comoving_nd.f90`