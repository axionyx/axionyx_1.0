
# Namespace eos\_module


[**Class List**](annotated.md) **>** [**eos\_module**](namespaceeos__module.md)


















## Public Attributes

| Type | Name |
| ---: | :--- |
|  logical, parameter, public | [**eos\_assume\_neutral**](namespaceeos__module.md#variable-eos-assume-neutral)   = = .false.<br> |
|  real(c\_double), public | [**vode\_atol\_scaled**](namespaceeos__module.md#variable-vode-atol-scaled)  <br> |
|  real(c\_double), public | [**vode\_rtol**](namespaceeos__module.md#variable-vode-rtol)  <br> |
|  real(rt), public | [**xacc**](namespaceeos__module.md#variable-xacc)  <br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**eos**](namespaceeos__module.md#function-eos) (input input, dens dens, temp temp, xmass xmass, pres pres, eint eint, c\_v c\_v, dPdT dPdT, dPdR dPdR, dEdT dEdT, entropy entropy, do\_eos\_diag do\_eos\_diag) <br> |
|  subroutine, public | [**eos\_init\_small\_pres**](namespaceeos__module.md#function-eos-init-small-pres) (R R, T T, Ne Ne, P P, a a) <br> |
|  subroutine | [**fort\_setup\_eos\_params**](namespaceeos__module.md#function-fort-setup-eos-params) (xacc\_in xacc\_in, vode\_rtol\_in vode\_rtol\_in, vode\_atol\_scaled\_in vode\_atol\_scaled\_in) <br> |
|  subroutine | [**ion\_n\_vec**](namespaceeos__module.md#function-ion-n-vec) (JH JH, JHe JHe, U U, nh nh, ne ne, nhp nhp, nhep nhep, nhepp nhepp, t t, vec\_count vec\_count) <br> |
|  subroutine, public | [**iterate\_ne**](namespaceeos__module.md#function-iterate-ne) (JH JH, JHe JHe, z z, U U, t t, nh nh, ne ne, nh0 nh0, nhp nhp, nhe0 nhe0, nhep nhep, nhepp nhepp) <br> |
|  subroutine, public | [**iterate\_ne\_vec**](namespaceeos__module.md#function-iterate-ne-vec) (z z, U U, t t, nh nh, ne ne, nh0 nh0, nhp nhp, nhe0 nhe0, nhep nhep, nhepp nhepp, veclen veclen) <br> |
|  subroutine, public | [**nyx\_eos\_given\_rt**](namespaceeos__module.md#function-nyx-eos-given-rt) (e e, P P, R R, T T, Ne Ne, a a) <br> |
|  subroutine, public | [**nyx\_eos\_given\_rt**](namespaceeos__module.md#function-nyx-eos-given-rt) (e e, P P, R R, T T, Ne Ne, comoving\_a comoving\_a) <br> |
|  subroutine, public | [**nyx\_eos\_given\_rt\_vec**](namespaceeos__module.md#function-nyx-eos-given-rt-vec) (e e, P P, R R, T T, Ne Ne, a a, veclen veclen) <br> |
|  subroutine, public | [**nyx\_eos\_nh0\_and\_nhep**](namespaceeos__module.md#function-nyx-eos-nh0-and-nhep) (JH JH, JHe JHe, z z, rho rho, e e, nh0 nh0, nhep nhep) <br> |
|  subroutine | [**nyx\_eos\_s\_given\_re**](namespaceeos__module.md#function-nyx-eos-s-given-re) (S S, R R, T T, Ne Ne, a a) <br> |
|  subroutine, public | [**nyx\_eos\_s\_given\_re**](namespaceeos__module.md#function-nyx-eos-s-given-re) (S S, R R, e e, T T, Ne Ne, comoving\_a comoving\_a) <br> |
|  subroutine | [**nyx\_eos\_soundspeed**](namespaceeos__module.md#function-nyx-eos-soundspeed) (c c, R R, e e) <br> |
|  subroutine, public | [**nyx\_eos\_t\_given\_re**](namespaceeos__module.md#function-nyx-eos-t-given-re) (JH JH, JHe JHe, T T, Ne Ne, R\_in R\_in, e\_in e\_in, a a, species species) <br> |
|  subroutine, public | [**nyx\_eos\_t\_given\_re**](namespaceeos__module.md#function-nyx-eos-t-given-re) (JH JH, JHe JHe, T T, Ne Ne, R R, e e, comoving\_a comoving\_a) <br> |
|  subroutine, public | [**nyx\_eos\_t\_given\_re\_vec**](namespaceeos__module.md#function-nyx-eos-t-given-re-vec) (T T, Ne Ne, R\_in R\_in, e\_in e\_in, a a, veclen veclen) <br> |








## Public Attributes Documentation


### <a href="#variable-eos-assume-neutral" id="variable-eos-assume-neutral">variable eos\_assume\_neutral </a>


```cpp
logical, parameter, public eos_module::eos_assume_neutral;
```



### <a href="#variable-vode-atol-scaled" id="variable-vode-atol-scaled">variable vode\_atol\_scaled </a>


```cpp
real(c_double), public eos_module::vode_atol_scaled;
```



### <a href="#variable-vode-rtol" id="variable-vode-rtol">variable vode\_rtol </a>


```cpp
real(c_double), public eos_module::vode_rtol;
```



### <a href="#variable-xacc" id="variable-xacc">variable xacc </a>


```cpp
real(rt), public eos_module::xacc;
```


## Public Functions Documentation


### <a href="#function-eos" id="function-eos">function eos </a>


```cpp
subroutine, public eos_module::eos (
    input input,
    dens dens,
    temp temp,
    xmass xmass,
    pres pres,
    eint eint,
    c_v c_v,
    dPdT dPdT,
    dPdR dPdR,
    dEdT dEdT,
    entropy entropy,
    do_eos_diag do_eos_diag
) 
```



### <a href="#function-eos-init-small-pres" id="function-eos-init-small-pres">function eos\_init\_small\_pres </a>


```cpp
subroutine, public eos_module::eos_init_small_pres (
    R R,
    T T,
    Ne Ne,
    P P,
    a a
) 
```



### <a href="#function-fort-setup-eos-params" id="function-fort-setup-eos-params">function fort\_setup\_eos\_params </a>


```cpp
subroutine eos_module::fort_setup_eos_params (
    xacc_in xacc_in,
    vode_rtol_in vode_rtol_in,
    vode_atol_scaled_in vode_atol_scaled_in
) 
```



### <a href="#function-ion-n-vec" id="function-ion-n-vec">function ion\_n\_vec </a>


```cpp
subroutine eos_module::ion_n_vec (
    JH JH,
    JHe JHe,
    U U,
    nh nh,
    ne ne,
    nhp nhp,
    nhep nhep,
    nhepp nhepp,
    t t,
    vec_count vec_count
) 
```



### <a href="#function-iterate-ne" id="function-iterate-ne">function iterate\_ne </a>


```cpp
subroutine, public eos_module::iterate_ne (
    JH JH,
    JHe JHe,
    z z,
    U U,
    t t,
    nh nh,
    ne ne,
    nh0 nh0,
    nhp nhp,
    nhe0 nhe0,
    nhep nhep,
    nhepp nhepp
) 
```



### <a href="#function-iterate-ne-vec" id="function-iterate-ne-vec">function iterate\_ne\_vec </a>


```cpp
subroutine, public eos_module::iterate_ne_vec (
    z z,
    U U,
    t t,
    nh nh,
    ne ne,
    nh0 nh0,
    nhp nhp,
    nhe0 nhe0,
    nhep nhep,
    nhepp nhepp,
    veclen veclen
) 
```



### <a href="#function-nyx-eos-given-rt" id="function-nyx-eos-given-rt">function nyx\_eos\_given\_rt </a>


```cpp
subroutine, public eos_module::nyx_eos_given_rt (
    e e,
    P P,
    R R,
    T T,
    Ne Ne,
    a a
) 
```



### <a href="#function-nyx-eos-given-rt" id="function-nyx-eos-given-rt">function nyx\_eos\_given\_rt </a>


```cpp
subroutine, public eos_module::nyx_eos_given_rt (
    e e,
    P P,
    R R,
    T T,
    Ne Ne,
    comoving_a comoving_a
) 
```



### <a href="#function-nyx-eos-given-rt-vec" id="function-nyx-eos-given-rt-vec">function nyx\_eos\_given\_rt\_vec </a>


```cpp
subroutine, public eos_module::nyx_eos_given_rt_vec (
    e e,
    P P,
    R R,
    T T,
    Ne Ne,
    a a,
    veclen veclen
) 
```



### <a href="#function-nyx-eos-nh0-and-nhep" id="function-nyx-eos-nh0-and-nhep">function nyx\_eos\_nh0\_and\_nhep </a>


```cpp
subroutine, public eos_module::nyx_eos_nh0_and_nhep (
    JH JH,
    JHe JHe,
    z z,
    rho rho,
    e e,
    nh0 nh0,
    nhep nhep
) 
```



### <a href="#function-nyx-eos-s-given-re" id="function-nyx-eos-s-given-re">function nyx\_eos\_s\_given\_re </a>


```cpp
subroutine eos_module::nyx_eos_s_given_re (
    S S,
    R R,
    T T,
    Ne Ne,
    a a
) 
```



### <a href="#function-nyx-eos-s-given-re" id="function-nyx-eos-s-given-re">function nyx\_eos\_s\_given\_re </a>


```cpp
subroutine, public eos_module::nyx_eos_s_given_re (
    S S,
    R R,
    e e,
    T T,
    Ne Ne,
    comoving_a comoving_a
) 
```



### <a href="#function-nyx-eos-soundspeed" id="function-nyx-eos-soundspeed">function nyx\_eos\_soundspeed </a>


```cpp
subroutine eos_module::nyx_eos_soundspeed (
    c c,
    R R,
    e e
) 
```



### <a href="#function-nyx-eos-t-given-re" id="function-nyx-eos-t-given-re">function nyx\_eos\_t\_given\_re </a>


```cpp
subroutine, public eos_module::nyx_eos_t_given_re (
    JH JH,
    JHe JHe,
    T T,
    Ne Ne,
    R_in R_in,
    e_in e_in,
    a a,
    species species
) 
```



### <a href="#function-nyx-eos-t-given-re" id="function-nyx-eos-t-given-re">function nyx\_eos\_t\_given\_re </a>


```cpp
subroutine, public eos_module::nyx_eos_t_given_re (
    JH JH,
    JHe JHe,
    T T,
    Ne Ne,
    R R,
    e e,
    comoving_a comoving_a
) 
```



### <a href="#function-nyx-eos-t-given-re-vec" id="function-nyx-eos-t-given-re-vec">function nyx\_eos\_t\_given\_re\_vec </a>


```cpp
subroutine, public eos_module::nyx_eos_t_given_re_vec (
    T T,
    Ne Ne,
    R_in R_in,
    e_in e_in,
    a a,
    veclen veclen
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/EOS/eos_hc.f90`