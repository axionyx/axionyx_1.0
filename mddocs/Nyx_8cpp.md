
# File Nyx.cpp


[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Nyx.cpp**](Nyx_8cpp.md)

[Go to the source code of this file.](Nyx_8cpp_source.md)



* `#include <iomanip>`
* `#include <algorithm>`
* `#include <vector>`
* `#include <iostream>`
* `#include <string>`
* `#include <unistd.h>`
* `#include <AMReX_CONSTANTS.H>`
* `#include <Nyx.H>`
* `#include <Nyx_F.H>`
* `#include <Derive_F.H>`
* `#include <AMReX_VisMF.H>`
* `#include <AMReX_TagBox.H>`
* `#include <AMReX_Utility.H>`
* `#include <AMReX_Print.H>`













## Public Attributes

| Type | Name |
| ---: | :--- |
|  const int | [**GimletSignal**](Nyx_8cpp.md#variable-gimletsignal)   = = 55<br> |
|  const int | [**NyxHaloFinderSignal**](Nyx_8cpp.md#variable-nyxhalofindersignal)   = = 42<br> |
|  int | [**simd\_width**](Nyx_8cpp.md#variable-simd-width)   = = 1<br> |
|  std::string | [**slice\_file**](Nyx_8cpp.md#variable-slice-file)   = = "slice\_"<br> |

## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  Real | [**dt\_cutoff**](Nyx_8cpp.md#variable-dt-cutoff)   = =  0<br> |
|  Real | [**fixed\_dt**](Nyx_8cpp.md#variable-fixed-dt)   = = -1.0<br> |
|  Real | [**initial\_dt**](Nyx_8cpp.md#variable-initial-dt)   = = -1.0<br> |
|  int | [**slice\_int**](Nyx_8cpp.md#variable-slice-int)   = = -1<br> |
|  int | [**slice\_nfiles**](Nyx_8cpp.md#variable-slice-nfiles)   = = 128<br> |
|  int | [**sum\_interval**](Nyx_8cpp.md#variable-sum-interval)   = = -1<br> |

## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**get\_comp\_e\_int**](Nyx_8cpp.md#function-get-comp-e-int) () <br> |
|  int | [**get\_comp\_temp**](Nyx_8cpp.md#function-get-comp-temp) () <br> |
|  int | [**get\_comp\_urho**](Nyx_8cpp.md#function-get-comp-urho) () <br> |
|  int | [**gimlet\_int**](Nyx_8cpp.md#function-gimlet-int) (0) <br> |
|  int | [**reeber\_int**](Nyx_8cpp.md#function-reeber-int) (0) <br> |








## Public Attributes Documentation


### <a href="#variable-gimletsignal" id="variable-gimletsignal">variable GimletSignal </a>


```cpp
const int GimletSignal;
```



### <a href="#variable-nyxhalofindersignal" id="variable-nyxhalofindersignal">variable NyxHaloFinderSignal </a>


```cpp
const int NyxHaloFinderSignal;
```



### <a href="#variable-simd-width" id="variable-simd-width">variable simd\_width </a>


```cpp
int simd_width;
```



### <a href="#variable-slice-file" id="variable-slice-file">variable slice\_file </a>


```cpp
std::string slice_file;
```


## Public Static Attributes Documentation


### <a href="#variable-dt-cutoff" id="variable-dt-cutoff">variable dt\_cutoff </a>


```cpp
Real dt_cutoff;
```



### <a href="#variable-fixed-dt" id="variable-fixed-dt">variable fixed\_dt </a>


```cpp
Real fixed_dt;
```



### <a href="#variable-initial-dt" id="variable-initial-dt">variable initial\_dt </a>


```cpp
Real initial_dt;
```



### <a href="#variable-slice-int" id="variable-slice-int">variable slice\_int </a>


```cpp
int slice_int;
```



### <a href="#variable-slice-nfiles" id="variable-slice-nfiles">variable slice\_nfiles </a>


```cpp
int slice_nfiles;
```



### <a href="#variable-sum-interval" id="variable-sum-interval">variable sum\_interval </a>


```cpp
int sum_interval;
```


## Public Functions Documentation


### <a href="#function-get-comp-e-int" id="function-get-comp-e-int">function get\_comp\_e\_int </a>


```cpp
int get_comp_e_int () 
```



### <a href="#function-get-comp-temp" id="function-get-comp-temp">function get\_comp\_temp </a>


```cpp
int get_comp_temp () 
```



### <a href="#function-get-comp-urho" id="function-get-comp-urho">function get\_comp\_urho </a>


```cpp
int get_comp_urho () 
```



### <a href="#function-gimlet-int" id="function-gimlet-int">function gimlet\_int </a>


```cpp
int gimlet_int (
    0
) 
```



### <a href="#function-reeber-int" id="function-reeber-int">function reeber\_int </a>


```cpp
int reeber_int (
    0
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/Nyx.cpp`