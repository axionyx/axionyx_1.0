
# File AMReX\_Amr.cpp


[**File List**](files.md) **>** [**AMReX\_axionyx**](dir_5c77c3c750fcf9b051dca9dbb6924de0.md) **>** [**AMReX\_Amr.cpp**](AMReX__Amr_8cpp.md)

[Go to the source code of this file.](AMReX__Amr_8cpp_source.md)



* `#include <algorithm>`
* `#include <cstdio>`
* `#include <list>`
* `#include <iostream>`
* `#include <iomanip>`
* `#include <sstream>`
* `#include <limits>`
* `#include <cmath>`
* `#include <sys/types.h>`
* `#include <sys/stat.h>`
* `#include <unistd.h>`
* `#include <AMReX_Geometry.H>`
* `#include <AMReX_TagBox.H>`
* `#include <AMReX_Array.H>`
* `#include <AMReX_Vector.H>`
* `#include <AMReX_CoordSys.H>`
* `#include <AMReX_ParmParse.H>`
* `#include <AMReX_BoxDomain.H>`
* `#include <AMReX_Cluster.H>`
* `#include <AMReX_LevelBld.H>`
* `#include <AMReX_AmrLevel.H>`
* `#include <AMReX_PROB_AMR_F.H>`
* `#include <AMReX_Amr.H>`
* `#include <AMReX_ParallelDescriptor.H>`
* `#include <AMReX_Utility.H>`
* `#include <AMReX_DistributionMapping.H>`
* `#include <AMReX_FabSet.H>`
* `#include <AMReX_StateData.H>`
* `#include <AMReX_PlotFileUtil.H>`
* `#include <AMReX_Print.H>`









## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**amrex**](namespaceamrex.md) <br> |
| namespace | [**anonymous\_namespace{AMReX\_Amr.cpp}**](namespaceamrex_1_1anonymous__namespace_02AMReX__Amr_8cpp_03.md) <br> |













## Macros

| Type | Name |
| ---: | :--- |
| define  | [**STRIP**](AMReX__Amr_8cpp.md#define-strip)  () while( is.get() != '\n' ) {}<br> |
| define  | [**STRIP**](AMReX__Amr_8cpp.md#define-strip)  () while( is.get() != '\n' ) {}<br> |
## Macro Definition Documentation



### <a href="#define-strip" id="define-strip">define STRIP </a>


```cpp
#define STRIP () 
```



### <a href="#define-strip" id="define-strip">define STRIP </a>


```cpp
#define STRIP () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/AMReX_axionyx/AMReX_Amr.cpp`