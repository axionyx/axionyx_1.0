
# File Gravity.cpp


[**File List**](files.md) **>** [**Gravity**](dir_fdbf5007869eac89a42b1cd44aeda050.md) **>** [**Gravity.cpp**](Gravity_8cpp.md)

[Go to the source code of this file.](Gravity_8cpp_source.md)



* `#include <cmath>`
* `#include <AMReX_ParmParse.H>`
* `#include "Gravity.H"`
* `#include "Nyx.H"`
* `#include <Gravity_F.H>`
* `#include <Nyx_F.H>`
* `#include <AMReX_MultiGrid.H>`
* `#include <AMReX_Laplacian.H>`
* `#include <AMReX_MacBndry.H>`
* `#include <AMReX_LO_BCTYPES.H>`
* `#include <AMReX_MLMG.H>`
* `#include <AMReX_MLPoisson.H>`














## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  Real | [**Ggravity**](Gravity_8cpp.md#variable-ggravity)   = = 0<br> |

## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**fort\_get\_grav\_const**](Gravity_8cpp.md#function-fort-get-grav-const) (Real \* Gconst) <br> |







## Macros

| Type | Name |
| ---: | :--- |
| define  | [**MAX\_LEV**](Gravity_8cpp.md#define-max-lev)  () 15<br> |

## Public Static Attributes Documentation


### <a href="#variable-ggravity" id="variable-ggravity">variable Ggravity </a>


```cpp
Real Ggravity;
```


## Public Functions Documentation


### <a href="#function-fort-get-grav-const" id="function-fort-get-grav-const">function fort\_get\_grav\_const </a>


```cpp
void fort_get_grav_const (
    Real * Gconst
) 
```

## Macro Definition Documentation



### <a href="#define-max-lev" id="define-max-lev">define MAX\_LEV </a>


```cpp
#define MAX_LEV () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Gravity/Gravity.cpp`