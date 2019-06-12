
# File Nyx\_error\_F.H


[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Tagging**](dir_c14a965952b26c2f69053cc66c8fb69f.md) **>** [**Nyx\_error\_F.H**](Nyx__error__F_8H.md)

[Go to the source code of this file.](Nyx__error__F_8H_source.md)



* `#include <AMReX_BLFort.H>`













## Public Attributes

| Type | Name |
| ---: | :--- |
|  int const int const int const int const int const int const int const int const amrex::Real const int const amrex::Real \* | [**avg**](Nyx__error__F_8H.md#variable-avg)  <br> |
|  const int const int const int const int const amrex::Real const amrex::Real const amrex::Real const int const amrex::Real \* | [**avg\_den**](Nyx__error__F_8H.md#variable-avg-den)  <br> |
|  int const int const int \* | [**clearval**](Nyx__error__F_8H.md#variable-clearval)  <br> |
|  int const int const int const int const int const int const int const int | [**domhi**](Nyx__error__F_8H.md#variable-domhi)  <br> |
|  int const int const int const int const int const int const int | [**domlo**](Nyx__error__F_8H.md#variable-domlo)  <br> |
|  int const int const int const int const int const int const int const int const amrex::Real | [**dx**](Nyx__error__F_8H.md#variable-dx)  <br> |
|  int const int const int const int const int | [**hi**](Nyx__error__F_8H.md#variable-hi)  <br> |
|  int const int const int const int const int const int const int const int const amrex::Real const amrex::Real const amrex::Real const amrex::Real const int \* | [**level**](Nyx__error__F_8H.md#variable-level)  <br> |
|  int const int const int const int | [**lo**](Nyx__error__F_8H.md#variable-lo)  <br> |
|  int const int const int const int const int const int \* | [**ncomp**](Nyx__error__F_8H.md#variable-ncomp)  <br> |
|  int const int const int const int const int const int const int const int const amrex::Real const amrex::Real const amrex::Real | [**problo**](Nyx__error__F_8H.md#variable-problo)  <br> |
|  int \* | [**tag**](Nyx__error__F_8H.md#variable-tag)  <br> |
|  int const int \* | [**tagval**](Nyx__error__F_8H.md#variable-tagval)  <br> |
|  int const int const int const int const int const int const int const int const amrex::Real const amrex::Real const amrex::Real const amrex::Real \* | [**time**](Nyx__error__F_8H.md#variable-time)  <br> |
|  const int const int const int const int const amrex::Real const amrex::Real const amrex::Real const int const amrex::Real const int \* | [**trigger**](Nyx__error__F_8H.md#variable-trigger)  <br> |
|  int const int const int const int const int const int const int const int const amrex::Real const amrex::Real | [**xlo**](Nyx__error__F_8H.md#variable-xlo)  <br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**ARLIM\_P**](Nyx__error__F_8H.md#function-arlim-p) (tag\_lo) <br> |
|  int | [**ARLIM\_P**](Nyx__error__F_8H.md#function-arlim-p) (tag\_hi) <br> |
|  int const int const int | [**BL\_FORT\_FAB\_ARG**](Nyx__error__F_8H.md#function-bl-fort-fab-arg) (var) <br> |
|  int const int const int | [**BL\_FORT\_FAB\_ARG**](Nyx__error__F_8H.md#function-bl-fort-fab-arg) (den) <br> |
|  int const int const int | [**BL\_FORT\_FAB\_ARG**](Nyx__error__F_8H.md#function-bl-fort-fab-arg) (vel) <br> |
|  int const int const int | [**BL\_FORT\_FAB\_ARG**](Nyx__error__F_8H.md#function-bl-fort-fab-arg) (temp) <br> |
|  int const int const int | [**BL\_FORT\_FAB\_ARG**](Nyx__error__F_8H.md#function-bl-fort-fab-arg) (press) <br> |
|  int const int const int | [**BL\_FORT\_FAB\_ARG**](Nyx__error__F_8H.md#function-bl-fort-fab-arg) (ls) <br> |
|  int const int const int | [**BL\_FORT\_FAB\_ARG**](Nyx__error__F_8H.md#function-bl-fort-fab-arg) (state) <br> |
|  int const int const int | [**BL\_FORT\_FAB\_ARG**](Nyx__error__F_8H.md#function-bl-fort-fab-arg) (axions) <br> |








## Public Attributes Documentation


### <a href="#variable-avg" id="variable-avg">variable avg </a>


```cpp
int const int const int const int const int const int const int const int const amrex::Real const int const amrex::Real * avg;
```



### <a href="#variable-avg-den" id="variable-avg-den">variable avg\_den </a>


```cpp
const int const int const int const int const amrex::Real const amrex::Real const amrex::Real const int const amrex::Real* avg_den;
```



### <a href="#variable-clearval" id="variable-clearval">variable clearval </a>


```cpp
int const int const int * clearval;
```



### <a href="#variable-domhi" id="variable-domhi">variable domhi </a>


```cpp
int const int const int const int const int const int const int const int domhi;
```



### <a href="#variable-domlo" id="variable-domlo">variable domlo </a>


```cpp
int const int const int const int const int const int const int domlo;
```



### <a href="#variable-dx" id="variable-dx">variable dx </a>


```cpp
int const int const int const int const int const int const int const int const amrex::Real dx[];
```



### <a href="#variable-hi" id="variable-hi">variable hi </a>


```cpp
int const int const int const int const int hi[];
```



### <a href="#variable-level" id="variable-level">variable level </a>


```cpp
int const int const int const int const int const int const int const int const amrex::Real const int* level;
```



### <a href="#variable-lo" id="variable-lo">variable lo </a>


```cpp
int const int const int const int lo[];
```



### <a href="#variable-ncomp" id="variable-ncomp">variable ncomp </a>


```cpp
int const int const int const int const int const int* ncomp;
```



### <a href="#variable-problo" id="variable-problo">variable problo </a>


```cpp
const int const int const int const int const amrex::Real const amrex::Real const amrex::Real problo;
```



### <a href="#variable-tag" id="variable-tag">variable tag </a>


```cpp
int * tag;
```



### <a href="#variable-tagval" id="variable-tagval">variable tagval </a>


```cpp
int const int * tagval;
```



### <a href="#variable-time" id="variable-time">variable time </a>


```cpp
int const int const int const int const int const int const int const int const amrex::Real const amrex::Real const amrex::Real const amrex::Real* time;
```



### <a href="#variable-trigger" id="variable-trigger">variable trigger </a>


```cpp
const int const int const int const int const amrex::Real const amrex::Real const amrex::Real const int const amrex::Real const int* trigger;
```



### <a href="#variable-xlo" id="variable-xlo">variable xlo </a>


```cpp
const int const int const int const int const amrex::Real const amrex::Real xlo[];
```


## Public Functions Documentation


### <a href="#function-arlim-p" id="function-arlim-p">function ARLIM\_P </a>


```cpp
int ARLIM_P (
    tag_lo
) 
```



### <a href="#function-arlim-p" id="function-arlim-p">function ARLIM\_P </a>


```cpp
int ARLIM_P (
    tag_hi
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
int const int const int BL_FORT_FAB_ARG (
    var
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
int const int const int BL_FORT_FAB_ARG (
    den
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
int const int const int BL_FORT_FAB_ARG (
    vel
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
int const int const int BL_FORT_FAB_ARG (
    temp
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
int const int const int BL_FORT_FAB_ARG (
    press
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
int const int const int BL_FORT_FAB_ARG (
    ls
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
int const int const int BL_FORT_FAB_ARG (
    state
) 
```



### <a href="#function-bl-fort-fab-arg" id="function-bl-fort-fab-arg">function BL\_FORT\_FAB\_ARG </a>


```cpp
int const int const int BL_FORT_FAB_ARG (
    axions
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/Tagging/Nyx_error_F.H`