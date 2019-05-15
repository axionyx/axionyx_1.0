
# Namespace atomic\_rates\_module


[**Class List**](annotated.md) **>** [**atomic\_rates\_module**](namespaceatomic__rates__module.md)


















## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rt), dimension(ncooltab+1), public | [**alphad**](namespaceatomic__rates__module.md#variable-alphad)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**alphahep**](namespaceatomic__rates__module.md#variable-alphahep)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**alphahepp**](namespaceatomic__rates__module.md#variable-alphahepp)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**alphahp**](namespaceatomic__rates__module.md#variable-alphahp)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**betaff1**](namespaceatomic__rates__module.md#variable-betaff1)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**betaff4**](namespaceatomic__rates__module.md#variable-betaff4)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**betah0**](namespaceatomic__rates__module.md#variable-betah0)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**betahe0**](namespaceatomic__rates__module.md#variable-betahe0)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**betahep**](namespaceatomic__rates__module.md#variable-betahep)  <br> |
|  real(rt), parameter, public | [**boltzmann**](namespaceatomic__rates__module.md#variable-boltzmann)   = = 1.3806e-16<br> |
|  real(rt), parameter, public | [**deltat**](namespaceatomic__rates__module.md#variable-deltat)   = = (TCOOLMAX - TCOOLMIN)/NCOOLTAB<br> |
|  real(rt), save, public | [**eh0**](namespaceatomic__rates__module.md#variable-eh0)  <br> |
|  real(rt), save, public | [**ehe0**](namespaceatomic__rates__module.md#variable-ehe0)  <br> |
|  real(rt), save, public | [**ehep**](namespaceatomic__rates__module.md#variable-ehep)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**gammaeh0**](namespaceatomic__rates__module.md#variable-gammaeh0)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**gammaehe0**](namespaceatomic__rates__module.md#variable-gammaehe0)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**gammaehep**](namespaceatomic__rates__module.md#variable-gammaehep)  <br> |
|  real(rt), save, public | [**ggh0**](namespaceatomic__rates__module.md#variable-ggh0)  <br> |
|  real(rt), save, public | [**gghe0**](namespaceatomic__rates__module.md#variable-gghe0)  <br> |
|  real(rt), save, public | [**gghep**](namespaceatomic__rates__module.md#variable-gghep)  <br> |
|  real(rt), save, public | [**mean\_rhob**](namespaceatomic__rates__module.md#variable-mean-rhob)  <br> |
|  real(rt), parameter, public | [**mproton**](namespaceatomic__rates__module.md#variable-mproton)   = = 1.6726231d-24<br> |
|  integer, parameter, public | [**ncooltab**](namespaceatomic__rates__module.md#variable-ncooltab)   = =2000<br> |
|  real(rt), dimension(ncooltab+1), public | [**rechep**](namespaceatomic__rates__module.md#variable-rechep)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**rechepp**](namespaceatomic__rates__module.md#variable-rechepp)  <br> |
|  real(rt), dimension(ncooltab+1), public | [**rechp**](namespaceatomic__rates__module.md#variable-rechp)  <br> |
|  real(rt), parameter, public | [**tcoolmax**](namespaceatomic__rates__module.md#variable-tcoolmax)   = = 9.0d0<br> |
|  real(rt), parameter, public | [**tcoolmax\_r**](namespaceatomic__rates__module.md#variable-tcoolmax-r)   = = 10.0d0\*\*TCOOLMAX<br> |
|  real(rt), parameter, public | [**tcoolmin**](namespaceatomic__rates__module.md#variable-tcoolmin)   = = 0.0d0<br> |
|  real(rt), parameter, public | [**tcoolmin\_r**](namespaceatomic__rates__module.md#variable-tcoolmin-r)   = = 10.0d0\*\*TCOOLMIN<br> |
|  real(rt), save, public | [**this\_z**](namespaceatomic__rates__module.md#variable-this-z)  <br> |
|  real(rt), save, public | [**uvb\_density\_a**](namespaceatomic__rates__module.md#variable-uvb-density-a)   = = 1.0d0<br> |
|  real(rt), save, public | [**uvb\_density\_b**](namespaceatomic__rates__module.md#variable-uvb-density-b)   = = 0.0d0<br> |
|  real(rt), public | [**xhydrogen**](namespaceatomic__rates__module.md#variable-xhydrogen)   = = 0.76d0<br> |
|  real(rt), public | [**yhelium**](namespaceatomic__rates__module.md#variable-yhelium)   = = 7.8947368421d-2<br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**fort\_interp\_to\_this\_z**](namespaceatomic__rates__module.md#function-fort-interp-to-this-z) (z z) <br> |
|  subroutine | [**fort\_tabulate\_rates**](namespaceatomic__rates__module.md#function-fort-tabulate-rates) () <br> |








## Public Attributes Documentation


### <a href="#variable-alphad" id="variable-alphad">variable alphad </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::alphad;
```



### <a href="#variable-alphahep" id="variable-alphahep">variable alphahep </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::alphahep;
```



### <a href="#variable-alphahepp" id="variable-alphahepp">variable alphahepp </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::alphahepp;
```



### <a href="#variable-alphahp" id="variable-alphahp">variable alphahp </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::alphahp;
```



### <a href="#variable-betaff1" id="variable-betaff1">variable betaff1 </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::betaff1;
```



### <a href="#variable-betaff4" id="variable-betaff4">variable betaff4 </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::betaff4;
```



### <a href="#variable-betah0" id="variable-betah0">variable betah0 </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::betah0;
```



### <a href="#variable-betahe0" id="variable-betahe0">variable betahe0 </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::betahe0;
```



### <a href="#variable-betahep" id="variable-betahep">variable betahep </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::betahep;
```



### <a href="#variable-boltzmann" id="variable-boltzmann">variable boltzmann </a>


```cpp
real(rt), parameter, public atomic_rates_module::boltzmann;
```



### <a href="#variable-deltat" id="variable-deltat">variable deltat </a>


```cpp
real(rt), parameter, public atomic_rates_module::deltat;
```



### <a href="#variable-eh0" id="variable-eh0">variable eh0 </a>


```cpp
real(rt), save, public atomic_rates_module::eh0;
```



### <a href="#variable-ehe0" id="variable-ehe0">variable ehe0 </a>


```cpp
real(rt), save, public atomic_rates_module::ehe0;
```



### <a href="#variable-ehep" id="variable-ehep">variable ehep </a>


```cpp
real(rt), save, public atomic_rates_module::ehep;
```



### <a href="#variable-gammaeh0" id="variable-gammaeh0">variable gammaeh0 </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::gammaeh0;
```



### <a href="#variable-gammaehe0" id="variable-gammaehe0">variable gammaehe0 </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::gammaehe0;
```



### <a href="#variable-gammaehep" id="variable-gammaehep">variable gammaehep </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::gammaehep;
```



### <a href="#variable-ggh0" id="variable-ggh0">variable ggh0 </a>


```cpp
real(rt), save, public atomic_rates_module::ggh0;
```



### <a href="#variable-gghe0" id="variable-gghe0">variable gghe0 </a>


```cpp
real(rt), save, public atomic_rates_module::gghe0;
```



### <a href="#variable-gghep" id="variable-gghep">variable gghep </a>


```cpp
real(rt), save, public atomic_rates_module::gghep;
```



### <a href="#variable-mean-rhob" id="variable-mean-rhob">variable mean\_rhob </a>


```cpp
real(rt), save, public atomic_rates_module::mean_rhob;
```



### <a href="#variable-mproton" id="variable-mproton">variable mproton </a>


```cpp
real(rt), parameter, public atomic_rates_module::mproton;
```



### <a href="#variable-ncooltab" id="variable-ncooltab">variable ncooltab </a>


```cpp
integer, parameter, public atomic_rates_module::ncooltab;
```



### <a href="#variable-rechep" id="variable-rechep">variable rechep </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::rechep;
```



### <a href="#variable-rechepp" id="variable-rechepp">variable rechepp </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::rechepp;
```



### <a href="#variable-rechp" id="variable-rechp">variable rechp </a>


```cpp
real(rt), dimension(ncooltab+1), public atomic_rates_module::rechp;
```



### <a href="#variable-tcoolmax" id="variable-tcoolmax">variable tcoolmax </a>


```cpp
real(rt), parameter, public atomic_rates_module::tcoolmax;
```



### <a href="#variable-tcoolmax-r" id="variable-tcoolmax-r">variable tcoolmax\_r </a>


```cpp
real(rt), parameter, public atomic_rates_module::tcoolmax_r;
```



### <a href="#variable-tcoolmin" id="variable-tcoolmin">variable tcoolmin </a>


```cpp
real(rt), parameter, public atomic_rates_module::tcoolmin;
```



### <a href="#variable-tcoolmin-r" id="variable-tcoolmin-r">variable tcoolmin\_r </a>


```cpp
real(rt), parameter, public atomic_rates_module::tcoolmin_r;
```



### <a href="#variable-this-z" id="variable-this-z">variable this\_z </a>


```cpp
real(rt), save, public atomic_rates_module::this_z;
```



### <a href="#variable-uvb-density-a" id="variable-uvb-density-a">variable uvb\_density\_a </a>


```cpp
real(rt), save, public atomic_rates_module::uvb_density_a;
```



### <a href="#variable-uvb-density-b" id="variable-uvb-density-b">variable uvb\_density\_b </a>


```cpp
real(rt), save, public atomic_rates_module::uvb_density_b;
```



### <a href="#variable-xhydrogen" id="variable-xhydrogen">variable xhydrogen </a>


```cpp
real(rt), public atomic_rates_module::xhydrogen;
```



### <a href="#variable-yhelium" id="variable-yhelium">variable yhelium </a>


```cpp
real(rt), public atomic_rates_module::yhelium;
```


## Public Functions Documentation


### <a href="#function-fort-interp-to-this-z" id="function-fort-interp-to-this-z">function fort\_interp\_to\_this\_z </a>


```cpp
subroutine atomic_rates_module::fort_interp_to_this_z (
    z z
) 
```



### <a href="#function-fort-tabulate-rates" id="function-fort-tabulate-rates">function fort\_tabulate\_rates </a>


```cpp
subroutine atomic_rates_module::fort_tabulate_rates () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/EOS/atomic_rates.f90`