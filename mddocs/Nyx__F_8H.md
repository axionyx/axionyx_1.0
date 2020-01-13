
# File Nyx\_F.H


[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Nyx\_F.H**](Nyx__F_8H.md)

[Go to the source code of this file.](Nyx__F_8H_source.md)



* `#include <AMReX_BLFort.H>`
* `#include <AMReX_REAL.H>`















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**RhsFnReal**](Nyx__F_8H.md#function-rhsfnreal) (double t, double \* u, double \* udot, double \* rpar, int neq) <br> |
|  void | [**adjust\_heat\_cool**](Nyx__F_8H.md#function-adjust-heat-cool) (const int lo, const int hi, BL\_FORT\_FAB\_ARG(S\_old), BL\_FORT\_FAB\_ARG(S\_new), BL\_FORT\_FAB\_ARG(ext\_src\_old), BL\_FORT\_FAB\_ARG(ext\_src\_new), const amrex::Real \* a\_old, const amrex::Real \* a\_new, const amrex::Real \* dt) <br> |
|  void | [**denfill**](Nyx__F_8H.md#function-denfill) (BL\_FORT\_FAB\_ARG(state), const int dlo, const int dhi, const amrex::Real dx, const amrex::Real glo, const amrex::Real \* time, const int bc) <br> |
|  void | [**filcc**](Nyx__F_8H.md#function-filcc) (const amrex::Real \* q, ARLIM\_P(q\_lo), ARLIM\_P(q\_hi), const int \* domlo, const int \* domhi, const amrex::Real \* dx\_crse, const amrex::Real \* xlo, const int \* bc) <br> |
|  void | [**fort\_add\_grav\_source**](Nyx__F_8H.md#function-fort-add-grav-source) (const int lo, const int hi, const  BL\_FORT\_FAB\_ARG, BL\_FORT\_FAB\_ARG(u\_out), const  BL\_FORT\_FAB\_ARG, const amrex::Real \* dt, const amrex::Real \* a\_old, const amrex::Real \* a\_new) <br> |
|  void | [**fort\_alloc\_simd\_vec**](Nyx__F_8H.md#function-fort-alloc-simd-vec) () <br> |
|  void | [**fort\_avgdown**](Nyx__F_8H.md#function-fort-avgdown) (BL\_FORT\_FAB\_ARG(crse\_fab), const int & nc, const  BL\_FORT\_FAB\_ARG, const int ovlo, const int ovhi, const int rat) <br> |
|  void | [**fort\_check\_initial\_species**](Nyx__F_8H.md#function-fort-check-initial-species) (const int \* lo, const int \* hi, BL\_FORT\_FAB\_ARG(state)) <br> |
|  void | [**fort\_compute\_gas\_frac**](Nyx__F_8H.md#function-fort-compute-gas-frac) (const int lo, const int hi, const amrex::Real dx, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, amrex::Real \* rho\_ave, amrex::Real \* T\_cut, amrex::Real \* rho\_cut, amrex::Real \* whim\_mass, amrex::Real \* whim\_vol, amrex::Real \* hh\_mass, amrex::Real \* hh\_vol, amrex::Real \* igm\_mass, amrex::Real \* igm\_vol, amrex::Real \* mass\_sum, amrex::Real \* vol\_sum) <br> |
|  void | [**fort\_compute\_max\_temp\_loc**](Nyx__F_8H.md#function-fort-compute-max-temp-loc) (const int lo, const int hi, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const amrex::Real \* max\_temp, const amrex::Real \* den\_maxt, const int \* imax, const int \* jmax, const int \* kmax) <br> |
|  void | [**fort\_compute\_rho\_temp**](Nyx__F_8H.md#function-fort-compute-rho-temp) (const int lo, const int hi, const amrex::Real dx, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, amrex::Real \* rho\_ave, amrex::Real \* rho\_T\_sum, amrex::Real \* T\_sum, amrex::Real \* Tinv\_sum, amrex::Real \* T\_meanrho\_sum, amrex::Real \* rho\_sum, amrex::Real \* vol\_sum, amrex::Real \* vol\_mn\_sum) <br> |
|  void | [**fort\_compute\_temp**](Nyx__F_8H.md#function-fort-compute-temp) (const int lo, const int hi, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, amrex::Real \* comoving\_a, const int \* print\_fortran\_warnings) <br> |
|  void | [**fort\_compute\_temp\_vec**](Nyx__F_8H.md#function-fort-compute-temp-vec) (const int lo, const int hi, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, amrex::Real \* comoving\_a, const int \* print\_fortran\_warnings) <br> |
|  void | [**fort\_correct\_gsrc**](Nyx__F_8H.md#function-fort-correct-gsrc) (const int lo, const int hi, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, BL\_FORT\_FAB\_ARG(S\_new), const amrex::Real \* a\_old, const amrex::Real \* a\_new, const amrex::Real \* dt) <br> |
|  void | [**fort\_dealloc\_simd\_vec**](Nyx__F_8H.md#function-fort-dealloc-simd-vec) () <br> |
|  void | [**fort\_enforce\_consistent\_e**](Nyx__F_8H.md#function-fort-enforce-consistent-e) (const int \* lo, const int \* hi, BL\_FORT\_FAB\_ARG(state)) <br> |
|  void | [**fort\_enforce\_nonnegative\_species**](Nyx__F_8H.md#function-fort-enforce-nonnegative-species) (BL\_FORT\_FAB\_ARG(S\_new), const int lo, const int hi, const int \* print\_fortran\_warnings) <br> |
|  void | [**fort\_estdt**](Nyx__F_8H.md#function-fort-estdt) (const  BL\_FORT\_FAB\_ARG, const int lo, const int hi, const amrex::Real dx, amrex::Real \* dt, amrex::Real \* comoving\_a) <br> |
|  void | [**fort\_estdt\_comoving\_a**](Nyx__F_8H.md#function-fort-estdt-comoving-a) (amrex::Real \* old\_a, amrex::Real \* new\_dummy\_a, amrex::Real \* dt, amrex::Real \* change\_allowed, amrex::Real \* fixed\_da\_interval, amrex::Real \* final\_a, int \* dt\_modified) <br> |
|  void | [**fort\_ext\_src**](Nyx__F_8H.md#function-fort-ext-src) (const int \* lo, const int \* hi, BL\_FORT\_FAB\_ARG(old\_state), BL\_FORT\_FAB\_ARG(new\_state), BL\_FORT\_FAB\_ARG(old\_diag), BL\_FORT\_FAB\_ARG(new\_diag), BL\_FORT\_FAB\_ARG(ext\_src), const amrex::Real \* prob\_lo, const amrex::Real \* dx, const amrex::Real \* time, const amrex::Real \* z, const amrex::Real \* dt) <br> |
|  void | [**fort\_get\_aux\_names**](Nyx__F_8H.md#function-fort-get-aux-names) (int \* aux\_names, int \* iaux, int \* len) <br> |
|  void | [**fort\_get\_hubble**](Nyx__F_8H.md#function-fort-get-hubble) (amrex::Real \* hubble) <br> |
|  void | [**fort\_get\_method\_params**](Nyx__F_8H.md#function-fort-get-method-params) (int \* HYP\_GROW) <br> |
|  void | [**fort\_get\_num\_aux**](Nyx__F_8H.md#function-fort-get-num-aux) (int \* naux) <br> |
|  void | [**fort\_get\_num\_spec**](Nyx__F_8H.md#function-fort-get-num-spec) (int \* nspec) <br> |
|  void | [**fort\_get\_spec\_names**](Nyx__F_8H.md#function-fort-get-spec-names) (int \* spec\_names, int \* ispec, int \* len) <br> |
|  void | [**fort\_init\_e\_from\_rhoe**](Nyx__F_8H.md#function-fort-init-e-from-rhoe) (BL\_FORT\_FAB\_ARG(state), const int \* num\_state, const int \* lo, const int \* hi, amrex::Real \* comoving\_a) <br> |
|  void | [**fort\_init\_e\_from\_t**](Nyx__F_8H.md#function-fort-init-e-from-t) (BL\_FORT\_FAB\_ARG(state), const int \* num\_state, BL\_FORT\_FAB\_ARG(diag), const int \* num\_diag, const int \* lo, const int \* hi, amrex::Real \* comoving\_a) <br> |
|  void | [**fort\_init\_zhi**](Nyx__F_8H.md#function-fort-init-zhi) (const int \* lo, const int \* hi, const int & num\_diag, BL\_FORT\_FAB\_ARG(diag\_eos), const int & ratio, BL\_FORT\_FAB\_ARG(zhi)) <br> |
|  void | [**fort\_initdata**](Nyx__F_8H.md#function-fort-initdata) (const int & level, const amrex::Real & time, const int \* lo, const int \* hi, const int & num\_state, BL\_FORT\_FAB\_ARG(state), const int & num\_diag, BL\_FORT\_FAB\_ARG(diag\_eos), const amrex::Real dx, const amrex::Real xlo, const amrex::Real xhi, const int \* domlo, const int \* domhi) <br> |
|  void | [**fort\_integrate\_comoving\_a**](Nyx__F_8H.md#function-fort-integrate-comoving-a) (amrex::Real \* old\_a, amrex::Real \* new\_a, amrex::Real \* dt) <br> |
|  void | [**fort\_integrate\_comoving\_a\_to\_a**](Nyx__F_8H.md#function-fort-integrate-comoving-a-to-a) (amrex::Real \* old\_a, amrex::Real \* a\_value, amrex::Real \* dt) <br> |
|  void | [**fort\_integrate\_comoving\_a\_to\_z**](Nyx__F_8H.md#function-fort-integrate-comoving-a-to-z) (amrex::Real \* old\_a, amrex::Real \* z\_value, amrex::Real \* dt) <br> |
|  void | [**fort\_integrate\_time\_given\_a**](Nyx__F_8H.md#function-fort-integrate-time-given-a) (amrex::Real \* a0, amrex::Real \* a1, amrex::Real \* dt) <br> |
|  void | [**fort\_interp\_to\_this\_z**](Nyx__F_8H.md#function-fort-interp-to-this-z) (const amrex::Real \* z) <br> |
|  void | [**fort\_make\_hydro\_sources**](Nyx__F_8H.md#function-fort-make-hydro-sources) (const amrex::Real \* time, const int lo, const int hi, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, BL\_FORT\_FAB\_ARG(hydro\_src), BL\_FORT\_FAB\_ARG(divu\_cc), const  BL\_FORT\_FAB\_ARG, const amrex::Real dx, const amrex::Real \* dt, D\_DECL(const BL\_FORT\_FAB\_ARG(xflux), const BL\_FORT\_FAB\_ARG(yflux), const BL\_FORT\_FAB\_ARG(zflux)), const amrex::Real \* a\_old, const amrex::Real \* a\_new, const int \* print\_fortran\_warnings) <br> |
|  void | [**fort\_network\_init**](Nyx__F_8H.md#function-fort-network-init) () <br> |
|  void | [**fort\_ode\_eos\_finalize**](Nyx__F_8H.md#function-fort-ode-eos-finalize) (double \* e\_out, double \* rpar, int neq) <br> |
|  void | [**fort\_ode\_eos\_setup**](Nyx__F_8H.md#function-fort-ode-eos-setup) (const amrex::Real & a, const amrex::Real & half\_dt) <br> |
|  void | [**fort\_set\_eos\_params**](Nyx__F_8H.md#function-fort-set-eos-params) (const amrex::Real & h\_species\_in, const amrex::Real & he\_species\_in) <br> |
|  void | [**fort\_set\_hubble**](Nyx__F_8H.md#function-fort-set-hubble) (const amrex::Real & hubble) <br> |
|  void | [**fort\_set\_method\_params**](Nyx__F_8H.md#function-fort-set-method-params) (const int & dm, const int & NumAdv, const int & Ndiag, const int & do\_hydro, const int & ppm\_type, const int & ppm\_ref, const int & ppm\_flatten\_before\_integrals, const int & use\_colglaz, const int & use\_flattening, const int & corner\_coupling, const int & version\_2, const int & use\_const\_species, const amrex::Real & gamma\_in, const int & normalize\_species, const int & heat\_cool\_type, const int & inhomo\_reion, const int & use\_axions) <br> |
|  void | [**fort\_set\_omb**](Nyx__F_8H.md#function-fort-set-omb) (const amrex::Real & frac) <br> |
|  void | [**fort\_set\_omm**](Nyx__F_8H.md#function-fort-set-omm) (const amrex::Real & omm) <br> |
|  void | [**fort\_set\_problem\_params**](Nyx__F_8H.md#function-fort-set-problem-params) (const int & dm, const int \* physbc\_lo, const int \* physbc\_hi, const int & Outflow\_value, const int & Symmetry\_value, const int & coord\_type) <br> |
|  void | [**fort\_set\_small\_values**](Nyx__F_8H.md#function-fort-set-small-values) (const amrex::Real \* average\_dens, const amrex::Real \* average\_temp, const amrex::Real \* comoving\_a, const amrex::Real \* small\_dens, const amrex::Real \* small\_temp, const amrex::Real \* small\_pres) <br> |
|  void | [**fort\_set\_xhydrogen**](Nyx__F_8H.md#function-fort-set-xhydrogen) (amrex::Real & xhydrogen\_in) <br> |
|  void | [**fort\_setup\_eos\_params**](Nyx__F_8H.md#function-fort-setup-eos-params) (amrex::Real \* eos\_nr\_eps, amrex::Real \* vode\_rtol, amrex::Real \* vode\_atol\_scaled) <br> |
|  void | [**fort\_syncgsrc**](Nyx__F_8H.md#function-fort-syncgsrc) (const int lo, const int hi, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, BL\_FORT\_FAB\_ARG(src), const amrex::Real \* a\_new, const amrex::Real & dt) <br> |
|  void | [**fort\_tabulate\_rates**](Nyx__F_8H.md#function-fort-tabulate-rates) () <br> |
|  void | [**fort\_update\_eos**](Nyx__F_8H.md#function-fort-update-eos) (double dt, double \* u, double \* uout, double \* rpar) <br> |
|  void | [**fort\_update\_state**](Nyx__F_8H.md#function-fort-update-state) (const int lo, const int hi, const  BL\_FORT\_FAB\_ARG, BL\_FORT\_FAB\_ARG(u\_out), const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const amrex::Real \* dt, const amrex::Real \* a\_old, const amrex::Real \* a\_new, const int \* print\_fortran\_warnings) <br> |
|  void | [**generic\_fill**](Nyx__F_8H.md#function-generic-fill) (BL\_FORT\_FAB\_ARG(state), const int dlo, const int dhi, const amrex::Real dx, const amrex::Real glo, const amrex::Real \* time, const int bc) <br> |
|  void | [**get\_rhoe**](Nyx__F_8H.md#function-get-rhoe) (const int lo, const int hi, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG, const  BL\_FORT\_FAB\_ARG) <br> |
|  void | [**hypfill**](Nyx__F_8H.md#function-hypfill) (BL\_FORT\_FAB\_ARG(state), const int dlo, const int dhi, const amrex::Real dx, const amrex::Real glo, const amrex::Real \* time, const int bc) <br> |
|  void | [**integrate\_state**](Nyx__F_8H.md#function-integrate-state) (const int \* lo, const int \* hi, BL\_FORT\_FAB\_ARG(state), BL\_FORT\_FAB\_ARG(diag\_eos), const amrex::Real \* a, const amrex::Real \* delta\_time, const int \* min\_iter, const int \* max\_iter) <br> |
|  void | [**integrate\_state\_fcvode\_with\_source**](Nyx__F_8H.md#function-integrate-state-fcvode-with-source) (const int \* lo, const int \* hi, BL\_FORT\_FAB\_ARG(state\_old), BL\_FORT\_FAB\_ARG(state\_new), BL\_FORT\_FAB\_ARG(diag\_eos), BL\_FORT\_FAB\_ARG(hydro\_src), BL\_FORT\_FAB\_ARG(reset\_src), BL\_FORT\_FAB\_ARG(IR), const amrex::Real \* a, const amrex::Real \* delta\_time, const int \* min\_iter, const int \* max\_iter) <br> |
|  void | [**integrate\_state\_with\_source**](Nyx__F_8H.md#function-integrate-state-with-source) (const int \* lo, const int \* hi, BL\_FORT\_FAB\_ARG(state\_old), BL\_FORT\_FAB\_ARG(state\_new), BL\_FORT\_FAB\_ARG(diag\_eos), BL\_FORT\_FAB\_ARG(hydro\_src), BL\_FORT\_FAB\_ARG(reset\_src), BL\_FORT\_FAB\_ARG(IR), const amrex::Real \* a, const amrex::Real \* delta\_time, const int \* min\_iter, const int \* max\_iter) <br> |
|  void | [**reset\_internal\_e**](Nyx__F_8H.md#function-reset-internal-e) (const int lo, const int hi, BL\_FORT\_FAB\_ARG(S\_new), BL\_FORT\_FAB\_ARG(D\_new), BL\_FORT\_FAB\_ARG(reset\_e\_src), const int \* print\_fortran\_warnings, amrex::Real \* comoving\_a) <br> |
|  void | [**set\_simd**](Nyx__F_8H.md#function-set-simd) (const int \* simd\_width) <br> |
|  void | [**sum\_over\_level**](Nyx__F_8H.md#function-sum-over-level) (BL\_FORT\_FAB\_ARG(rho), const int lo, const int hi, const amrex::Real dx, amrex::Real \* sum) <br> |
|  void | [**sum\_prod\_prod**](Nyx__F_8H.md#function-sum-prod-prod) (BL\_FORT\_FAB\_ARG(fab1), BL\_FORT\_FAB\_ARG(fab2), BL\_FORT\_FAB\_ARG(fab3), const int lo, const int hi, const amrex::Real dx, amrex::Real \* sum) <br> |
|  void | [**sum\_product**](Nyx__F_8H.md#function-sum-product) (BL\_FORT\_FAB\_ARG(fab1), BL\_FORT\_FAB\_ARG(fab2), const int lo, const int hi, const amrex::Real dx, amrex::Real \* sum) <br> |
|  void | [**time\_center\_sources**](Nyx__F_8H.md#function-time-center-sources) (const int lo, const int hi, BL\_FORT\_FAB\_ARG(S\_new), BL\_FORT\_FAB\_ARG(ext\_src\_old), BL\_FORT\_FAB\_ARG(ext\_src\_new), const amrex::Real \* a\_old, const amrex::Real \* a\_new, const amrex::Real \* dt, const int \* print\_fortran\_warnings) <br> |








## Public Functions Documentation


### <a href="#function-rhsfnreal" id="function-rhsfnreal">function RhsFnReal </a>


```cpp
void RhsFnReal (
    double t,
    double * u,
    double * udot,
    double * rpar,
    int neq
) 
```



### <a href="#function-adjust-heat-cool" id="function-adjust-heat-cool">function adjust\_heat\_cool </a>


```cpp
void adjust_heat_cool (
    const int lo,
    const int hi,
    BL_FORT_FAB_ARG(S_old),
    BL_FORT_FAB_ARG(S_new),
    BL_FORT_FAB_ARG(ext_src_old),
    BL_FORT_FAB_ARG(ext_src_new),
    const amrex::Real * a_old,
    const amrex::Real * a_new,
    const amrex::Real * dt
) 
```



### <a href="#function-denfill" id="function-denfill">function denfill </a>


```cpp
void denfill (
    BL_FORT_FAB_ARG(state),
    const int dlo,
    const int dhi,
    const amrex::Real dx,
    const amrex::Real glo,
    const amrex::Real * time,
    const int bc
) 
```



### <a href="#function-filcc" id="function-filcc">function filcc </a>


```cpp
void filcc (
    const amrex::Real * q,
    ARLIM_P(q_lo),
    ARLIM_P(q_hi),
    const int * domlo,
    const int * domhi,
    const amrex::Real * dx_crse,
    const amrex::Real * xlo,
    const int * bc
) 
```



### <a href="#function-fort-add-grav-source" id="function-fort-add-grav-source">function fort\_add\_grav\_source </a>


```cpp
void fort_add_grav_source (
    const int lo,
    const int hi,
    const BL_FORT_FAB_ARG,
    BL_FORT_FAB_ARG(u_out),
    const BL_FORT_FAB_ARG,
    const amrex::Real * dt,
    const amrex::Real * a_old,
    const amrex::Real * a_new
) 
```



### <a href="#function-fort-alloc-simd-vec" id="function-fort-alloc-simd-vec">function fort\_alloc\_simd\_vec </a>


```cpp
void fort_alloc_simd_vec () 
```



### <a href="#function-fort-avgdown" id="function-fort-avgdown">function fort\_avgdown </a>


```cpp
void fort_avgdown (
    BL_FORT_FAB_ARG(crse_fab),
    const int & nc,
    const BL_FORT_FAB_ARG,
    const int ovlo,
    const int ovhi,
    const int rat
) 
```



### <a href="#function-fort-check-initial-species" id="function-fort-check-initial-species">function fort\_check\_initial\_species </a>


```cpp
void fort_check_initial_species (
    const int * lo,
    const int * hi,
    BL_FORT_FAB_ARG(state)
) 
```



### <a href="#function-fort-compute-gas-frac" id="function-fort-compute-gas-frac">function fort\_compute\_gas\_frac </a>


```cpp
void fort_compute_gas_frac (
    const int lo,
    const int hi,
    const amrex::Real dx,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    amrex::Real * rho_ave,
    amrex::Real * T_cut,
    amrex::Real * rho_cut,
    amrex::Real * whim_mass,
    amrex::Real * whim_vol,
    amrex::Real * hh_mass,
    amrex::Real * hh_vol,
    amrex::Real * igm_mass,
    amrex::Real * igm_vol,
    amrex::Real * mass_sum,
    amrex::Real * vol_sum
) 
```



### <a href="#function-fort-compute-max-temp-loc" id="function-fort-compute-max-temp-loc">function fort\_compute\_max\_temp\_loc </a>


```cpp
void fort_compute_max_temp_loc (
    const int lo,
    const int hi,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const amrex::Real * max_temp,
    const amrex::Real * den_maxt,
    const int * imax,
    const int * jmax,
    const int * kmax
) 
```



### <a href="#function-fort-compute-rho-temp" id="function-fort-compute-rho-temp">function fort\_compute\_rho\_temp </a>


```cpp
void fort_compute_rho_temp (
    const int lo,
    const int hi,
    const amrex::Real dx,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    amrex::Real * rho_ave,
    amrex::Real * rho_T_sum,
    amrex::Real * T_sum,
    amrex::Real * Tinv_sum,
    amrex::Real * T_meanrho_sum,
    amrex::Real * rho_sum,
    amrex::Real * vol_sum,
    amrex::Real * vol_mn_sum
) 
```



### <a href="#function-fort-compute-temp" id="function-fort-compute-temp">function fort\_compute\_temp </a>


```cpp
void fort_compute_temp (
    const int lo,
    const int hi,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    amrex::Real * comoving_a,
    const int * print_fortran_warnings
) 
```



### <a href="#function-fort-compute-temp-vec" id="function-fort-compute-temp-vec">function fort\_compute\_temp\_vec </a>


```cpp
void fort_compute_temp_vec (
    const int lo,
    const int hi,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    amrex::Real * comoving_a,
    const int * print_fortran_warnings
) 
```



### <a href="#function-fort-correct-gsrc" id="function-fort-correct-gsrc">function fort\_correct\_gsrc </a>


```cpp
void fort_correct_gsrc (
    const int lo,
    const int hi,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    BL_FORT_FAB_ARG(S_new),
    const amrex::Real * a_old,
    const amrex::Real * a_new,
    const amrex::Real * dt
) 
```



### <a href="#function-fort-dealloc-simd-vec" id="function-fort-dealloc-simd-vec">function fort\_dealloc\_simd\_vec </a>


```cpp
void fort_dealloc_simd_vec () 
```



### <a href="#function-fort-enforce-consistent-e" id="function-fort-enforce-consistent-e">function fort\_enforce\_consistent\_e </a>


```cpp
void fort_enforce_consistent_e (
    const int * lo,
    const int * hi,
    BL_FORT_FAB_ARG(state)
) 
```



### <a href="#function-fort-enforce-nonnegative-species" id="function-fort-enforce-nonnegative-species">function fort\_enforce\_nonnegative\_species </a>


```cpp
void fort_enforce_nonnegative_species (
    BL_FORT_FAB_ARG(S_new),
    const int lo,
    const int hi,
    const int * print_fortran_warnings
) 
```



### <a href="#function-fort-estdt" id="function-fort-estdt">function fort\_estdt </a>


```cpp
void fort_estdt (
    const BL_FORT_FAB_ARG,
    const int lo,
    const int hi,
    const amrex::Real dx,
    amrex::Real * dt,
    amrex::Real * comoving_a
) 
```



### <a href="#function-fort-estdt-comoving-a" id="function-fort-estdt-comoving-a">function fort\_estdt\_comoving\_a </a>


```cpp
void fort_estdt_comoving_a (
    amrex::Real * old_a,
    amrex::Real * new_dummy_a,
    amrex::Real * dt,
    amrex::Real * change_allowed,
    amrex::Real * fixed_da_interval,
    amrex::Real * final_a,
    int * dt_modified
) 
```



### <a href="#function-fort-ext-src" id="function-fort-ext-src">function fort\_ext\_src </a>


```cpp
void fort_ext_src (
    const int * lo,
    const int * hi,
    BL_FORT_FAB_ARG(old_state),
    BL_FORT_FAB_ARG(new_state),
    BL_FORT_FAB_ARG(old_diag),
    BL_FORT_FAB_ARG(new_diag),
    BL_FORT_FAB_ARG(ext_src),
    const amrex::Real * prob_lo,
    const amrex::Real * dx,
    const amrex::Real * time,
    const amrex::Real * z,
    const amrex::Real * dt
) 
```



### <a href="#function-fort-get-aux-names" id="function-fort-get-aux-names">function fort\_get\_aux\_names </a>


```cpp
void fort_get_aux_names (
    int * aux_names,
    int * iaux,
    int * len
) 
```



### <a href="#function-fort-get-hubble" id="function-fort-get-hubble">function fort\_get\_hubble </a>


```cpp
void fort_get_hubble (
    amrex::Real * hubble
) 
```



### <a href="#function-fort-get-method-params" id="function-fort-get-method-params">function fort\_get\_method\_params </a>


```cpp
void fort_get_method_params (
    int * HYP_GROW
) 
```



### <a href="#function-fort-get-num-aux" id="function-fort-get-num-aux">function fort\_get\_num\_aux </a>


```cpp
void fort_get_num_aux (
    int * naux
) 
```



### <a href="#function-fort-get-num-spec" id="function-fort-get-num-spec">function fort\_get\_num\_spec </a>


```cpp
void fort_get_num_spec (
    int * nspec
) 
```



### <a href="#function-fort-get-spec-names" id="function-fort-get-spec-names">function fort\_get\_spec\_names </a>


```cpp
void fort_get_spec_names (
    int * spec_names,
    int * ispec,
    int * len
) 
```



### <a href="#function-fort-init-e-from-rhoe" id="function-fort-init-e-from-rhoe">function fort\_init\_e\_from\_rhoe </a>


```cpp
void fort_init_e_from_rhoe (
    BL_FORT_FAB_ARG(state),
    const int * num_state,
    const int * lo,
    const int * hi,
    amrex::Real * comoving_a
) 
```



### <a href="#function-fort-init-e-from-t" id="function-fort-init-e-from-t">function fort\_init\_e\_from\_t </a>


```cpp
void fort_init_e_from_t (
    BL_FORT_FAB_ARG(state),
    const int * num_state,
    BL_FORT_FAB_ARG(diag),
    const int * num_diag,
    const int * lo,
    const int * hi,
    amrex::Real * comoving_a
) 
```



### <a href="#function-fort-init-zhi" id="function-fort-init-zhi">function fort\_init\_zhi </a>


```cpp
void fort_init_zhi (
    const int * lo,
    const int * hi,
    const int & num_diag,
    BL_FORT_FAB_ARG(diag_eos),
    const int & ratio,
    BL_FORT_FAB_ARG(zhi)
) 
```



### <a href="#function-fort-initdata" id="function-fort-initdata">function fort\_initdata </a>


```cpp
void fort_initdata (
    const int & level,
    const amrex::Real & time,
    const int * lo,
    const int * hi,
    const int & num_state,
    BL_FORT_FAB_ARG(state),
    const int & num_diag,
    BL_FORT_FAB_ARG(diag_eos),
    const amrex::Real dx,
    const amrex::Real xlo,
    const amrex::Real xhi,
    const int * domlo,
    const int * domhi
) 
```



### <a href="#function-fort-integrate-comoving-a" id="function-fort-integrate-comoving-a">function fort\_integrate\_comoving\_a </a>


```cpp
void fort_integrate_comoving_a (
    amrex::Real * old_a,
    amrex::Real * new_a,
    amrex::Real * dt
) 
```



### <a href="#function-fort-integrate-comoving-a-to-a" id="function-fort-integrate-comoving-a-to-a">function fort\_integrate\_comoving\_a\_to\_a </a>


```cpp
void fort_integrate_comoving_a_to_a (
    amrex::Real * old_a,
    amrex::Real * a_value,
    amrex::Real * dt
) 
```



### <a href="#function-fort-integrate-comoving-a-to-z" id="function-fort-integrate-comoving-a-to-z">function fort\_integrate\_comoving\_a\_to\_z </a>


```cpp
void fort_integrate_comoving_a_to_z (
    amrex::Real * old_a,
    amrex::Real * z_value,
    amrex::Real * dt
) 
```



### <a href="#function-fort-integrate-time-given-a" id="function-fort-integrate-time-given-a">function fort\_integrate\_time\_given\_a </a>


```cpp
void fort_integrate_time_given_a (
    amrex::Real * a0,
    amrex::Real * a1,
    amrex::Real * dt
) 
```



### <a href="#function-fort-interp-to-this-z" id="function-fort-interp-to-this-z">function fort\_interp\_to\_this\_z </a>


```cpp
void fort_interp_to_this_z (
    const amrex::Real * z
) 
```



### <a href="#function-fort-make-hydro-sources" id="function-fort-make-hydro-sources">function fort\_make\_hydro\_sources </a>


```cpp
void fort_make_hydro_sources (
    const amrex::Real * time,
    const int lo,
    const int hi,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    BL_FORT_FAB_ARG(hydro_src),
    BL_FORT_FAB_ARG(divu_cc),
    const BL_FORT_FAB_ARG,
    const amrex::Real dx,
    const amrex::Real * dt,
    D_DECL(const BL_FORT_FAB_ARG(xflux), const BL_FORT_FAB_ARG(yflux), const BL_FORT_FAB_ARG(zflux)),
    const amrex::Real * a_old,
    const amrex::Real * a_new,
    const int * print_fortran_warnings
) 
```



### <a href="#function-fort-network-init" id="function-fort-network-init">function fort\_network\_init </a>


```cpp
void fort_network_init () 
```



### <a href="#function-fort-ode-eos-finalize" id="function-fort-ode-eos-finalize">function fort\_ode\_eos\_finalize </a>


```cpp
void fort_ode_eos_finalize (
    double * e_out,
    double * rpar,
    int neq
) 
```



### <a href="#function-fort-ode-eos-setup" id="function-fort-ode-eos-setup">function fort\_ode\_eos\_setup </a>


```cpp
void fort_ode_eos_setup (
    const amrex::Real & a,
    const amrex::Real & half_dt
) 
```



### <a href="#function-fort-set-eos-params" id="function-fort-set-eos-params">function fort\_set\_eos\_params </a>


```cpp
void fort_set_eos_params (
    const amrex::Real & h_species_in,
    const amrex::Real & he_species_in
) 
```



### <a href="#function-fort-set-hubble" id="function-fort-set-hubble">function fort\_set\_hubble </a>


```cpp
void fort_set_hubble (
    const amrex::Real & hubble
) 
```



### <a href="#function-fort-set-method-params" id="function-fort-set-method-params">function fort\_set\_method\_params </a>


```cpp
void fort_set_method_params (
    const int & dm,
    const int & NumAdv,
    const int & Ndiag,
    const int & do_hydro,
    const int & ppm_type,
    const int & ppm_ref,
    const int & ppm_flatten_before_integrals,
    const int & use_colglaz,
    const int & use_flattening,
    const int & corner_coupling,
    const int & version_2,
    const int & use_const_species,
    const amrex::Real & gamma_in,
    const int & normalize_species,
    const int & heat_cool_type,
    const int & inhomo_reion,
    const int & use_axions
) 
```



### <a href="#function-fort-set-omb" id="function-fort-set-omb">function fort\_set\_omb </a>


```cpp
void fort_set_omb (
    const amrex::Real & frac
) 
```



### <a href="#function-fort-set-omm" id="function-fort-set-omm">function fort\_set\_omm </a>


```cpp
void fort_set_omm (
    const amrex::Real & omm
) 
```



### <a href="#function-fort-set-problem-params" id="function-fort-set-problem-params">function fort\_set\_problem\_params </a>


```cpp
void fort_set_problem_params (
    const int & dm,
    const int * physbc_lo,
    const int * physbc_hi,
    const int & Outflow_value,
    const int & Symmetry_value,
    const int & coord_type
) 
```



### <a href="#function-fort-set-small-values" id="function-fort-set-small-values">function fort\_set\_small\_values </a>


```cpp
void fort_set_small_values (
    const amrex::Real * average_dens,
    const amrex::Real * average_temp,
    const amrex::Real * comoving_a,
    const amrex::Real * small_dens,
    const amrex::Real * small_temp,
    const amrex::Real * small_pres
) 
```



### <a href="#function-fort-set-xhydrogen" id="function-fort-set-xhydrogen">function fort\_set\_xhydrogen </a>


```cpp
void fort_set_xhydrogen (
    amrex::Real & xhydrogen_in
) 
```



### <a href="#function-fort-setup-eos-params" id="function-fort-setup-eos-params">function fort\_setup\_eos\_params </a>


```cpp
void fort_setup_eos_params (
    amrex::Real * eos_nr_eps,
    amrex::Real * vode_rtol,
    amrex::Real * vode_atol_scaled
) 
```



### <a href="#function-fort-syncgsrc" id="function-fort-syncgsrc">function fort\_syncgsrc </a>


```cpp
void fort_syncgsrc (
    const int lo,
    const int hi,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    BL_FORT_FAB_ARG(src),
    const amrex::Real * a_new,
    const amrex::Real & dt
) 
```



### <a href="#function-fort-tabulate-rates" id="function-fort-tabulate-rates">function fort\_tabulate\_rates </a>


```cpp
void fort_tabulate_rates () 
```



### <a href="#function-fort-update-eos" id="function-fort-update-eos">function fort\_update\_eos </a>


```cpp
void fort_update_eos (
    double dt,
    double * u,
    double * uout,
    double * rpar
) 
```



### <a href="#function-fort-update-state" id="function-fort-update-state">function fort\_update\_state </a>


```cpp
void fort_update_state (
    const int lo,
    const int hi,
    const BL_FORT_FAB_ARG,
    BL_FORT_FAB_ARG(u_out),
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const amrex::Real * dt,
    const amrex::Real * a_old,
    const amrex::Real * a_new,
    const int * print_fortran_warnings
) 
```



### <a href="#function-generic-fill" id="function-generic-fill">function generic\_fill </a>


```cpp
void generic_fill (
    BL_FORT_FAB_ARG(state),
    const int dlo,
    const int dhi,
    const amrex::Real dx,
    const amrex::Real glo,
    const amrex::Real * time,
    const int bc
) 
```



### <a href="#function-get-rhoe" id="function-get-rhoe">function get\_rhoe </a>


```cpp
void get_rhoe (
    const int lo,
    const int hi,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG,
    const BL_FORT_FAB_ARG
) 
```



### <a href="#function-hypfill" id="function-hypfill">function hypfill </a>


```cpp
void hypfill (
    BL_FORT_FAB_ARG(state),
    const int dlo,
    const int dhi,
    const amrex::Real dx,
    const amrex::Real glo,
    const amrex::Real * time,
    const int bc
) 
```



### <a href="#function-integrate-state" id="function-integrate-state">function integrate\_state </a>


```cpp
void integrate_state (
    const int * lo,
    const int * hi,
    BL_FORT_FAB_ARG(state),
    BL_FORT_FAB_ARG(diag_eos),
    const amrex::Real * a,
    const amrex::Real * delta_time,
    const int * min_iter,
    const int * max_iter
) 
```



### <a href="#function-integrate-state-fcvode-with-source" id="function-integrate-state-fcvode-with-source">function integrate\_state\_fcvode\_with\_source </a>


```cpp
void integrate_state_fcvode_with_source (
    const int * lo,
    const int * hi,
    BL_FORT_FAB_ARG(state_old),
    BL_FORT_FAB_ARG(state_new),
    BL_FORT_FAB_ARG(diag_eos),
    BL_FORT_FAB_ARG(hydro_src),
    BL_FORT_FAB_ARG(reset_src),
    BL_FORT_FAB_ARG(IR),
    const amrex::Real * a,
    const amrex::Real * delta_time,
    const int * min_iter,
    const int * max_iter
) 
```



### <a href="#function-integrate-state-with-source" id="function-integrate-state-with-source">function integrate\_state\_with\_source </a>


```cpp
void integrate_state_with_source (
    const int * lo,
    const int * hi,
    BL_FORT_FAB_ARG(state_old),
    BL_FORT_FAB_ARG(state_new),
    BL_FORT_FAB_ARG(diag_eos),
    BL_FORT_FAB_ARG(hydro_src),
    BL_FORT_FAB_ARG(reset_src),
    BL_FORT_FAB_ARG(IR),
    const amrex::Real * a,
    const amrex::Real * delta_time,
    const int * min_iter,
    const int * max_iter
) 
```



### <a href="#function-reset-internal-e" id="function-reset-internal-e">function reset\_internal\_e </a>


```cpp
void reset_internal_e (
    const int lo,
    const int hi,
    BL_FORT_FAB_ARG(S_new),
    BL_FORT_FAB_ARG(D_new),
    BL_FORT_FAB_ARG(reset_e_src),
    const int * print_fortran_warnings,
    amrex::Real * comoving_a
) 
```



### <a href="#function-set-simd" id="function-set-simd">function set\_simd </a>


```cpp
void set_simd (
    const int * simd_width
) 
```



### <a href="#function-sum-over-level" id="function-sum-over-level">function sum\_over\_level </a>


```cpp
void sum_over_level (
    BL_FORT_FAB_ARG(rho),
    const int lo,
    const int hi,
    const amrex::Real dx,
    amrex::Real * sum
) 
```



### <a href="#function-sum-prod-prod" id="function-sum-prod-prod">function sum\_prod\_prod </a>


```cpp
void sum_prod_prod (
    BL_FORT_FAB_ARG(fab1),
    BL_FORT_FAB_ARG(fab2),
    BL_FORT_FAB_ARG(fab3),
    const int lo,
    const int hi,
    const amrex::Real dx,
    amrex::Real * sum
) 
```



### <a href="#function-sum-product" id="function-sum-product">function sum\_product </a>


```cpp
void sum_product (
    BL_FORT_FAB_ARG(fab1),
    BL_FORT_FAB_ARG(fab2),
    const int lo,
    const int hi,
    const amrex::Real dx,
    amrex::Real * sum
) 
```



### <a href="#function-time-center-sources" id="function-time-center-sources">function time\_center\_sources </a>


```cpp
void time_center_sources (
    const int lo,
    const int hi,
    BL_FORT_FAB_ARG(S_new),
    BL_FORT_FAB_ARG(ext_src_old),
    BL_FORT_FAB_ARG(ext_src_new),
    const amrex::Real * a_old,
    const amrex::Real * a_new,
    const amrex::Real * dt,
    const int * print_fortran_warnings
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Nyx_F.H`