
# File Nyx.H


[**File List**](files.md) **>** [**Exec**](dir_43a12cefb7942b6f49b5b628aafd3192.md) **>** [**Test\_Only\_Axions**](dir_eb24725df855cf6c732a19e4912f662a.md) **>** [**Nyx.H**](Exec_2Test__Only__Axions_2Nyx_8H.md)

[Go to the source code of this file.](Exec_2Test__Only__Axions_2Nyx_8H_source.md)



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
| class | [**Nyx**](classNyx.md) <br>_AmrLevel-derived class for hyperbolic conservation equations for stellar media._  |

## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**StateType**](Exec_2Test__Only__Axions_2Nyx_8H.md#enum-statetype)  <br> |
| typedef [**NyxParticleContainer**](classNyxParticleContainer.md)&lt; 1+BL\_SPACEDIM+3 &gt; | [**StellarParticleContainer**](Exec_2Test__Only__Axions_2Nyx_8H.md#typedef-stellarparticlecontainer)  <br> |


## Public Attributes

| Type | Name |
| ---: | :--- |
|  int | [**gimlet\_int**](Exec_2Test__Only__Axions_2Nyx_8H.md#variable-gimlet-int)  <br>_time step interval for doing Gimlet post-processing_  |
|  int | [**reeber\_int**](Exec_2Test__Only__Axions_2Nyx_8H.md#variable-reeber-int)  <br>_time step interval for finding halos_  |










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
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Exec/Test_Only_Axions/Nyx.H`