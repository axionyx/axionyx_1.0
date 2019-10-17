
# File integrate\_state\_box\_3d.cpp


[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**integrate\_state\_box\_3d.cpp**](integrate__state__box__3d_8cpp.md)

[Go to the source code of this file.](integrate__state__box__3d_8cpp_source.md)



* `#include <fstream>`
* `#include <iomanip>`
* `#include <AMReX_ParmParse.H>`
* `#include <AMReX_Geometry.H>`
* `#include <AMReX_MultiFab.H>`
* `#include <AMReX_Print.H>`
* `#include <AMReX_PlotFileUtil.H>`
* `#include <AMReX_BLFort.H>`
* `#include <Nyx.H>`
* `#include <Nyx_F.H>`
* `#include <cvode/cvode.h>`
* `#include <cvode/cvode_diag.h>`
* `#include <sundials/sundials_types.h>`
* `#include <sundials/sundials_math.h>`
* `#include <nvector/nvector_serial.h>`
















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  void | [**PrintFinalStats**](integrate__state__box__3d_8cpp.md#function-printfinalstats) (void \* cvode\_mem) <br> |
|  void | [**PrintOutput**](integrate__state__box__3d_8cpp.md#function-printoutput) (realtype t, realtype umax, long int nst) <br> |
|  int | [**check\_retval**](integrate__state__box__3d_8cpp.md#function-check-retval) (void \* flagvalue, const char \* funcname, int opt) <br> |
|  int | [**f**](integrate__state__box__3d_8cpp.md#function-f) (realtype t, N\_Vector u, N\_Vector udot, void \* user\_data) <br> |







## Public Static Functions Documentation


### <a href="#function-printfinalstats" id="function-printfinalstats">function PrintFinalStats </a>


```cpp
static void PrintFinalStats (
    void * cvode_mem
) 
```



### <a href="#function-printoutput" id="function-printoutput">function PrintOutput </a>


```cpp
static void PrintOutput (
    realtype t,
    realtype umax,
    long int nst
) 
```



### <a href="#function-check-retval" id="function-check-retval">function check\_retval </a>


```cpp
static int check_retval (
    void * flagvalue,
    const char * funcname,
    int opt
) 
```



### <a href="#function-f" id="function-f">function f </a>


```cpp
static int f (
    realtype t,
    N_Vector u,
    N_Vector udot,
    void * user_data
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/HeatCool/integrate_state_box_3d.cpp`