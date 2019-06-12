
# File nyx\_main.cpp


[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**nyx\_main.cpp**](nyx__main_8cpp.md)

[Go to the source code of this file.](nyx__main_8cpp_source.md)



* `#include <iostream>`
* `#include <iomanip>`
* `#include <sstream>`
* `#include <unistd.h>`
* `#include <AMReX_CArena.H>`
* `#include <AMReX_REAL.H>`
* `#include <AMReX_Utility.H>`
* `#include <AMReX_IntVect.H>`
* `#include <AMReX_Box.H>`
* `#include <AMReX_Amr.H>`
* `#include <AMReX_ParmParse.H>`
* `#include <AMReX_ParallelDescriptor.H>`
* `#include <AMReX_AmrLevel.H>`
* `#include <AMReX_Geometry.H>`
* `#include <AMReX_MultiFab.H>`
* `#include <Nyx.H>`
* `#include "Nyx_output.H"`













## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::string | [**inputs\_name**](nyx__main_8cpp.md#variable-inputs-name)   = = ""<br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|  const int | [**GimletSignal**](nyx__main_8cpp.md#function-gimletsignal) (55) <br> |
|  const int | [**NyxHaloFinderSignal**](nyx__main_8cpp.md#function-nyxhalofindersignal) (42) <br> |
|  void | [**nyx\_main**](nyx__main_8cpp.md#function-nyx-main) (int argc, char \* argv) <br> |
|  const int | [**quitSignal**](nyx__main_8cpp.md#function-quitsignal) (- 44) <br> |
|  const int | [**resizeSignal**](nyx__main_8cpp.md#function-resizesignal) (43) <br> |








## Public Attributes Documentation


### <a href="#variable-inputs-name" id="variable-inputs-name">variable inputs\_name </a>


```cpp
std::string inputs_name;
```


## Public Functions Documentation


### <a href="#function-gimletsignal" id="function-gimletsignal">function GimletSignal </a>


```cpp
const int GimletSignal (
    55
) 
```



### <a href="#function-nyxhalofindersignal" id="function-nyxhalofindersignal">function NyxHaloFinderSignal </a>


```cpp
const int NyxHaloFinderSignal (
    42
) 
```



### <a href="#function-nyx-main" id="function-nyx-main">function nyx\_main </a>


```cpp
void nyx_main (
    int argc,
    char * argv
) 
```



### <a href="#function-quitsignal" id="function-quitsignal">function quitSignal </a>


```cpp
const int quitSignal (
    - 44
) 
```



### <a href="#function-resizesignal" id="function-resizesignal">function resizeSignal </a>


```cpp
const int resizeSignal (
    43
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/nyx_main.cpp`