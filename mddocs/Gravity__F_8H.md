
# File Gravity\_F.H


[**File List**](files.md) **>** [**Gravity**](dir_fdbf5007869eac89a42b1cd44aeda050.md) **>** [**Gravity\_F.H**](Gravity__F_8H.md)

[Go to the source code of this file.](Gravity__F_8H_source.md)



* `#include <AMReX_BLFort.H>`













## Public Attributes

| Type | Name |
| ---: | :--- |
|  const int const int const int const int const int \* | [**dir**](Gravity__F_8H.md#variable-dir)  <br> |
|  const int const int const int const int \* | [**domain\_hi**](Gravity__F_8H.md#variable-domain-hi)  <br> |
|  const int const int const int \* | [**domain\_lo**](Gravity__F_8H.md#variable-domain-lo)  <br> |
|  const int const int const const const const amrex::Real \* | [**dx**](Gravity__F_8H.md#variable-dx)  <br> |
|  const int const int \* | [**fhi**](Gravity__F_8H.md#variable-fhi)  <br> |
|  const int \* | [**flo**](Gravity__F_8H.md#variable-flo)  <br> |
|  const int const int \* | [**hi**](Gravity__F_8H.md#variable-hi)  <br> |
|  const int \* | [**lo**](Gravity__F_8H.md#variable-lo)  <br> |
|  const int const int const int \* | [**nc**](Gravity__F_8H.md#variable-nc)  <br> |
|  const int const int const int const int const int \* | [**nparticles**](Gravity__F_8H.md#variable-nparticles)  <br> |
|  const const int const int | [**ovhi**](Gravity__F_8H.md#variable-ovhi)  <br> |
|  const const int | [**ovlo**](Gravity__F_8H.md#variable-ovlo)  <br> |
|  const int const int const int const int const int const amrex::Real \* | [**part\_locs**](Gravity__F_8H.md#variable-part-locs)  <br> |
|  const int const int const int const int const int const amrex::Real const amrex::Real \* | [**part\_mass**](Gravity__F_8H.md#variable-part-mass)  <br> |
|  const const int const int const int | [**rat**](Gravity__F_8H.md#variable-rat)  <br> |
|  const int const int const int const int \* | [**refRatio**](Gravity__F_8H.md#variable-refratio)  <br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BL\_FORT\_FAB\_ARG**](Gravity__F_8H.md#function-bl-fort-fab-arg) (crse\_fab) <br> |
|  const | [**BL\_FORT\_FAB\_ARG**](Gravity__F_8H.md#function-bl-fort-fab-arg) (fine\_fab) <br> |
|  const int const int const int const int const int | [**BL\_FORT\_FAB\_ARG**](Gravity__F_8H.md#function-bl-fort-fab-arg) (fine) <br> |
|  const int const int const int const int const int const | [**BL\_FORT\_FAB\_ARG**](Gravity__F_8H.md#function-bl-fort-fab-arg) (crse) <br> |
|  const int const int | [**BL\_FORT\_FAB\_ARG**](Gravity__F_8H.md#function-bl-fort-fab-arg) (phi\_cc) <br> |
|  const int const int const | [**BL\_FORT\_FAB\_ARG**](Gravity__F_8H.md#function-bl-fort-fab-arg) (xgrad) <br> |
|  const int const int const const | [**BL\_FORT\_FAB\_ARG**](Gravity__F_8H.md#function-bl-fort-fab-arg) (ygrad) <br> |
|  const int const int const const const | [**BL\_FORT\_FAB\_ARG**](Gravity__F_8H.md#function-bl-fort-fab-arg) (zgrad) <br> |
|  const int const int const int const int | [**BL\_FORT\_FAB\_ARG**](Gravity__F_8H.md#function-bl-fort-fab-arg) (phi) <br> |








## Public Attributes Documentation


### <a href="#variable-dir" id="variable-dir">variable dir </a>


```cpp
const int const int const int const int const int * dir;
```



### <a href="#variable-domain-hi" id="variable-domain-hi">variable domain\_hi </a>


```cpp
const int const int const int const int* domain_hi;
```



### <a href="#variable-domain-lo" id="variable-domain-lo">variable domain\_lo </a>


```cpp
const int const int const int* domain_lo;
```



### <a href="#variable-dx" id="variable-dx">variable dx </a>


```cpp
int const int const int const int const int const int const int const int const amrex::Real dx;
```



### <a href="#variable-fhi" id="variable-fhi">variable fhi </a>


```cpp
const int const int* fhi;
```



### <a href="#variable-flo" id="variable-flo">variable flo </a>


```cpp
const int* flo;
```



### <a href="#variable-hi" id="variable-hi">variable hi </a>


```cpp
const int const int* hi;
```



### <a href="#variable-lo" id="variable-lo">variable lo </a>


```cpp
const int* lo;
```



### <a href="#variable-nc" id="variable-nc">variable nc </a>


```cpp
const int const int const int * nc;
```



### <a href="#variable-nparticles" id="variable-nparticles">variable nparticles </a>


```cpp
const int const int const int const int const int* nparticles;
```



### <a href="#variable-ovhi" id="variable-ovhi">variable ovhi </a>


```cpp
const const int const int ovhi[];
```



### <a href="#variable-ovlo" id="variable-ovlo">variable ovlo </a>


```cpp
const const int ovlo[];
```



### <a href="#variable-part-locs" id="variable-part-locs">variable part\_locs </a>


```cpp
const int const int const int const int const int const amrex::Real* part_locs;
```



### <a href="#variable-part-mass" id="variable-part-mass">variable part\_mass </a>


```cpp
const int const int const int const int const int const amrex::Real const amrex::Real* part_mass;
```



### <a href="#variable-rat" id="variable-rat">variable rat </a>


```cpp
const const int const int const int rat[];
```



### <a href="#variable-refratio" id="variable-refratio">variable refRatio </a>


```cpp
const int const int const int const int * refRatio;
```


## Public Functions Documentation


### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
BL_FORT_FAB_ARG (
    crse_fab
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
const BL_FORT_FAB_ARG (
    fine_fab
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
const int const int const int const int const int BL_FORT_FAB_ARG (
    fine
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
const int const int const int const int const int const BL_FORT_FAB_ARG (
    crse
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
const int const int BL_FORT_FAB_ARG (
    phi_cc
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
const int const int const BL_FORT_FAB_ARG (
    xgrad
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
const int const int const const BL_FORT_FAB_ARG (
    ygrad
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
const int const int const const const BL_FORT_FAB_ARG (
    zgrad
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
const int const int const int const int BL_FORT_FAB_ARG (
    phi
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Gravity/Gravity_F.H`