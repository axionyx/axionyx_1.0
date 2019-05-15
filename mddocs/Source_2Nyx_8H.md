
# File Nyx.H


[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Nyx.H**](Source_2Nyx_8H.md)

[Go to the source code of this file.](Source_2Nyx_8H_source.md)



* `#include <AMReX_BC_TYPES.H>`
* `#include <AMReX_AmrLevel.H>`
* `#include <AMReX_ErrorList.H>`
* `#include <AMReX_FluxRegister.H>`
* `#include "NyxParticleContainer.H"`
* `#include "DarkMatterParticleContainer.H"`
* `#include <iostream>`










## Classes

| Type | Name |
| ---: | :--- |
| class | [**Nyx**](classNyx.md) <br> |

## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**StateType**](Source_2Nyx_8H.md#enum-statetype)  <br> |
| typedef [**NyxParticleContainer**](classNyxParticleContainer.md)&lt; 1+BL\_SPACEDIM+3 &gt; | [**StellarParticleContainer**](Source_2Nyx_8H.md#typedef-stellarparticlecontainer)  <br> |


## Public Attributes

| Type | Name |
| ---: | :--- |
|  int | [**gimlet\_int**](Source_2Nyx_8H.md#variable-gimlet-int)  <br> |
|  int | [**reeber\_int**](Source_2Nyx_8H.md#variable-reeber-int)  <br> |










## Public Types Documentation


### <a href="#enum-statetype" id="enum-statetype">enum StateType </a>


```cpp
enum StateType {
    State_Type = 0,
    DiagEOS_Type,
    NUM_STATE_TYPE,
    State_Type = 0,
    DiagEOS_Type,
    NUM_STATE_TYPE
};
```



### <a href="#typedef-stellarparticlecontainer" id="typedef-stellarparticlecontainer">typedef StellarParticleContainer </a>


```cpp
typedef NyxParticleContainer<1+BL_SPACEDIM+3> StellarParticleContainer;
```


## Public Attributes Documentation


### <a href="#variable-gimlet-int" id="variable-gimlet-int">variable gimlet\_int </a>


```cpp
int gimlet_int;
```



### <a href="#variable-reeber-int" id="variable-reeber-int">variable reeber\_int </a>


```cpp
int reeber_int;
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/Nyx.H`