
# File Forcing.H


[**File List**](files.md) **>** [**Forcing**](dir_45682215f16eaf57f766b3c547de68bc.md) **>** [**Forcing.H**](Forcing_8H.md)

[Go to the source code of this file.](Forcing_8H_source.md)



* `#include <fstream>`
* `#include <AMReX_AmrLevel.H>`










## Classes

| Type | Name |
| ---: | :--- |
| class | [**StochasticForcing**](classStochasticForcing.md) <br> |

## Public Types

| Type | Name |
| ---: | :--- |
| typedef int | [**spect\_profile\_type**](Forcing_8H.md#typedef-spect-profile-type)  <br> |


## Public Attributes

| Type | Name |
| ---: | :--- |
|  const spect\_profile\_type | [**Band**](Forcing_8H.md#variable-band)   = = 2<br> |
|  const spect\_profile\_type | [**None**](Forcing_8H.md#variable-none)   = = 0<br> |
|  const spect\_profile\_type | [**Parabolic**](Forcing_8H.md#variable-parabolic)   = = 3<br> |
|  const spect\_profile\_type | [**Plane**](Forcing_8H.md#variable-plane)   = = 1<br> |









## Macros

| Type | Name |
| ---: | :--- |
| define  | [**MAX\_DIMENSION**](Forcing_8H.md#define-max-dimension)  () 3<br> |

## Public Types Documentation


### <a href="#typedef-spect-profile-type" id="typedef-spect-profile-type">typedef spect\_profile\_type </a>


```cpp
typedef int spect_profile_type;
```


## Public Attributes Documentation


### <a href="#variable-band" id="variable-band">variable Band </a>


```cpp
const spect_profile_type Band;
```



### <a href="#variable-none" id="variable-none">variable None </a>


```cpp
const spect_profile_type None;
```



### <a href="#variable-parabolic" id="variable-parabolic">variable Parabolic </a>


```cpp
const spect_profile_type Parabolic;
```



### <a href="#variable-plane" id="variable-plane">variable Plane </a>


```cpp
const spect_profile_type Plane;
```

## Macro Definition Documentation



### <a href="#define-max-dimension" id="define-max-dimension">define MAX\_DIMENSION </a>


```cpp
#define MAX_DIMENSION () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Forcing/Forcing.H`