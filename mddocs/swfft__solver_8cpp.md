
# File swfft\_solver.cpp


[**File List**](files.md) **>** [**Exec**](dir_43a12cefb7942b6f49b5b628aafd3192.md) **>** [**Test\_Only\_Axions**](dir_eb24725df855cf6c732a19e4912f662a.md) **>** [**swfft\_solver.cpp**](swfft__solver_8cpp.md)

[Go to the source code of this file.](swfft__solver_8cpp_source.md)



* `#include <AMReX_MultiFabUtil.H>`
* `#include <AMReX_VisMF.H>`
* `#include <AMReX_ParmParse.H>`
* `#include <Distribution.H>`
* `#include <AlignedAllocator.h>`
* `#include <Dfft.H>`
* `#include <string>`















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**swfft\_solver**](swfft__solver_8cpp.md#function-swfft-solver) (MultiFab & rhs, MultiFab & soln, Geometry & geom, int verbose) <br> |







## Macros

| Type | Name |
| ---: | :--- |
| define  | [**ALIGN**](swfft__solver_8cpp.md#define-align)  () 16<br> |

## Public Functions Documentation


### <a href="#function-swfft-solver" id="function-swfft-solver">function swfft\_solver </a>


```cpp
void swfft_solver (
    MultiFab & rhs,
    MultiFab & soln,
    Geometry & geom,
    int verbose
) 
```

## Macro Definition Documentation



### <a href="#define-align" id="define-align">define ALIGN </a>


```cpp
#define ALIGN () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Exec/Test_Only_Axions/swfft_solver.cpp`