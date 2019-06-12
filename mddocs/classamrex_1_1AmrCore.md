
# Class amrex::AmrCore


[**Class List**](annotated.md) **>** [**amrex**](namespaceamrex.md) **>** [**AmrCore**](classamrex_1_1AmrCore.md)



_Provide basic functionalities to set up an AMR hierarchy._ [More...](#detailed-description)

* `#include <AMReX_AmrCore.H>`



Inherits the following classes: AmrMesh


Inherited by the following classes: [amrex::Amr](classamrex_1_1Amr.md)










## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**AmrCore**](classamrex_1_1AmrCore.md#function-amrcore-1-3) () <br> |
|   | [**AmrCore**](classamrex_1_1AmrCore.md#function-amrcore-2-3) (const RealBox \* rb, int max\_level\_in, const Vector&lt; int &gt; & n\_cell\_in, int coord=-1, Vector&lt; IntVect &gt; ref\_ratios=Vector&lt; IntVect &gt;()) <br> |
|   | [**AmrCore**](classamrex_1_1AmrCore.md#function-amrcore-3-3) (const [**AmrCore**](classamrex_1_1AmrCore.md) & rhs) = delete<br> |
|  void | [**InitFromScratch**](classamrex_1_1AmrCore.md#function-initfromscratch) (Real time) <br>_Initialize BoxArray, DistributionMapping and data from scratch. Calling this function requires the derive class implement its own MakeNewLevelFromScratch to allocate and intialize data. Also note usually one needs to average the fine data down to coarse level after this._  |
|  int | [**Verbose**](classamrex_1_1AmrCore.md#function-verbose) () noexcept const<br> |
|  [**AmrCore**](classamrex_1_1AmrCore.md) & | [**operator=**](classamrex_1_1AmrCore.md#function-operator) (const [**AmrCore**](classamrex_1_1AmrCore.md) & rhs) = delete<br> |
|  void | [**printGridSummary**](classamrex_1_1AmrCore.md#function-printgridsummary) (std::ostream & os, int min\_lev, int max\_lev) noexcept const<br> |
| virtual void | [**regrid**](classamrex_1_1AmrCore.md#function-regrid) (int lbase, int iteration, Real time, bool initial=false) <br>_Rebuild levels finer than lbase._  |
| virtual  | [**~AmrCore**](classamrex_1_1AmrCore.md#function-amrcore) () <br> |

## Public Static Functions

| Type | Name |
| ---: | :--- |
|  void | [**Finalize**](classamrex_1_1AmrCore.md#function-finalize) () <br> |
|  void | [**Initialize**](classamrex_1_1AmrCore.md#function-initialize) () <br> |



## Protected Attributes

| Type | Name |
| ---: | :--- |
|  int | [**verbose**](classamrex_1_1AmrCore.md#variable-verbose)  <br> |


## Protected Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**ClearLevel**](classamrex_1_1AmrCore.md#function-clearlevel) (int lev) = 0<br>_Delete level data._  |
| virtual void | [**ErrorEst**](classamrex_1_1AmrCore.md#function-errorest) (int lev, TagBoxArray & tags, Real time, int ngrow) override = 0<br>_Tag cells for refinement. TagBoxArray tags is built on level lev grids._  |
| virtual void | [**MakeNewLevelFromCoarse**](classamrex_1_1AmrCore.md#function-makenewlevelfromcoarse) (int lev, Real time, const BoxArray & ba, const DistributionMapping & dm) = 0<br>_Make a new level using provided BoxArray and DistributionMapping and fill with interpolated coarse level data._  |
| virtual void | [**MakeNewLevelFromScratch**](classamrex_1_1AmrCore.md#function-makenewlevelfromscratch) (int lev, Real time, const BoxArray & ba, const DistributionMapping & dm) override = 0<br> |
| virtual void | [**RemakeLevel**](classamrex_1_1AmrCore.md#function-remakelevel) (int lev, Real time, const BoxArray & ba, const DistributionMapping & dm) = 0<br>_Remake an existing level using provided BoxArray and DistributionMapping and fill with existing fine and coarse data._  |


## Public Functions Documentation


### <a href="#function-amrcore-1-3" id="function-amrcore-1-3">function AmrCore [1/3]</a>


```cpp
amrex::AmrCore::AmrCore () 
```



### <a href="#function-amrcore-2-3" id="function-amrcore-2-3">function AmrCore [2/3]</a>


```cpp
amrex::AmrCore::AmrCore (
    const RealBox * rb,
    int max_level_in,
    const Vector< int > & n_cell_in,
    int coord=-1,
    Vector< IntVect > ref_ratios=Vector< IntVect >()
) 
```



### <a href="#function-amrcore-3-3" id="function-amrcore-3-3">function AmrCore [3/3]</a>


```cpp
amrex::AmrCore::AmrCore (
    const AmrCore & rhs
) = delete
```



### <a href="#function-initfromscratch" id="function-initfromscratch">function InitFromScratch </a>


```cpp
void amrex::AmrCore::InitFromScratch (
    Real time
) 
```



### <a href="#function-verbose" id="function-verbose">function Verbose </a>


```cpp
inline int amrex::AmrCore::Verbose () noexcept const
```



### <a href="#function-operator" id="function-operator">function operator= </a>


```cpp
AmrCore & amrex::AmrCore::operator= (
    const AmrCore & rhs
) = delete
```



### <a href="#function-printgridsummary" id="function-printgridsummary">function printGridSummary </a>


```cpp
void amrex::AmrCore::printGridSummary (
    std::ostream & os,
    int min_lev,
    int max_lev
) noexcept const
```



### <a href="#function-regrid" id="function-regrid">function regrid </a>


```cpp
virtual void amrex::AmrCore::regrid (
    int lbase,
    int iteration,
    Real time,
    bool initial=false
) 
```



### <a href="#function-amrcore" id="function-amrcore">function ~AmrCore </a>


```cpp
virtual amrex::AmrCore::~AmrCore () 
```


## Public Static Functions Documentation


### <a href="#function-finalize" id="function-finalize">function Finalize </a>


```cpp
static void amrex::AmrCore::Finalize () 
```



### <a href="#function-initialize" id="function-initialize">function Initialize </a>


```cpp
static void amrex::AmrCore::Initialize () 
```


## Protected Attributes Documentation


### <a href="#variable-verbose" id="variable-verbose">variable verbose </a>


```cpp
int amrex::AmrCore::verbose;
```


## Protected Functions Documentation


### <a href="#function-clearlevel" id="function-clearlevel">function ClearLevel </a>


```cpp
virtual void amrex::AmrCore::ClearLevel (
    int lev
) = 0
```



### <a href="#function-errorest" id="function-errorest">function ErrorEst </a>


```cpp
virtual void amrex::AmrCore::ErrorEst (
    int lev,
    TagBoxArray & tags,
    Real time,
    int ngrow
) override = 0
```



### <a href="#function-makenewlevelfromcoarse" id="function-makenewlevelfromcoarse">function MakeNewLevelFromCoarse </a>


```cpp
virtual void amrex::AmrCore::MakeNewLevelFromCoarse (
    int lev,
    Real time,
    const BoxArray & ba,
    const DistributionMapping & dm
) = 0
```



### <a href="#function-makenewlevelfromscratch" id="function-makenewlevelfromscratch">function MakeNewLevelFromScratch </a>


```cpp
virtual void amrex::AmrCore::MakeNewLevelFromScratch (
    int lev,
    Real time,
    const BoxArray & ba,
    const DistributionMapping & dm
) override = 0
```


Make a new level from scratch using provided BoxArray and DistributionMapping. Only used during initialization. 


        

### <a href="#function-remakelevel" id="function-remakelevel">function RemakeLevel </a>


```cpp
virtual void amrex::AmrCore::RemakeLevel (
    int lev,
    Real time,
    const BoxArray & ba,
    const DistributionMapping & dm
) = 0
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/AMReX_axionyx/AMReX_AmrCore.H`