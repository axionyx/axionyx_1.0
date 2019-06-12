
# Class amrex::AmrLevel


[**Class List**](annotated.md) **>** [**amrex**](namespaceamrex.md) **>** [**AmrLevel**](classamrex_1_1AmrLevel.md)



_Virtual base class for managing individual levels._ [_**AmrLevel**_](classamrex_1_1AmrLevel.md) _functions both as a container for state data on a level and also manages the advancement of data in time._

* `#include <AMReX_AmrLevel.H>`





Inherited by the following classes: [Nyx](classNyx.md),  [Nyx](classNyx.md)






## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**TimeLevel**](classamrex_1_1AmrLevel.md#enum-timelevel)  <br>_What time are we at?_  |




## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**CreateLevelDirectory**](classamrex_1_1AmrLevel.md#function-createleveldirectory) (const std::string & dir) <br>_Create the Level\_ directory for checkpoint and plot files._  |
|  const DistributionMapping & | [**DistributionMap**](classamrex_1_1AmrLevel.md#function-distributionmap) () noexcept const<br> |
|  const Box & | [**Domain**](classamrex_1_1AmrLevel.md#function-domain) () noexcept const<br>_Returns the indices defining physical domain._  |
|  const FabFactory&lt; FArrayBox &gt; & | [**Factory**](classamrex_1_1AmrLevel.md#function-factory) () noexcept const<br> |
|  void | [**FillCoarsePatch**](classamrex_1_1AmrLevel.md#function-fillcoarsepatch) (MultiFab & dest, int dcomp, Real time, int state\_idx, int scomp, int ncomp, int nghost=0) <br>_Interpolate from coarse level to the valid area in dest._  |
|  const Geometry & | [**Geom**](classamrex_1_1AmrLevel.md#function-geom) () noexcept const<br>_Returns the geometry object._  |
|  int | [**Level**](classamrex_1_1AmrLevel.md#function-level) () noexcept const<br>_Returns this_ [_**AmrLevel**_](classamrex_1_1AmrLevel.md) _._ |
|  void | [**LevelDirectoryNames**](classamrex_1_1AmrLevel.md#function-leveldirectorynames) (const std::string & dir, std::string & LevelDir, std::string & FullPath) <br>_Get the level directory names._  |
|  void | [**SetLevelDirectoryCreated**](classamrex_1_1AmrLevel.md#function-setleveldirectorycreated) (bool ldc) noexcept<br>_Set if the Level\_ directory was created or to clear the value. CreateLevelDirectory sets levelDirectoryCreated = true._  |
|  void | [**UpdateDistributionMaps**](classamrex_1_1AmrLevel.md#function-updatedistributionmaps) (DistributionMapping & dmap) <br>_Update the distribution maps in StateData based on the size of the map._  |
| virtual int | [**WorkEstType**](classamrex_1_1AmrLevel.md#function-workesttype) () <br>_Which state data type is for work estimates? -1 means none._  |
| virtual Real | [**advance**](classamrex_1_1AmrLevel.md#function-advance) (Real time, Real dt, int iteration, int ncycle) = 0<br>_Do an integration step on this level. Returns maximum safe time step. This is a pure virtual function and hence MUST be implemented by derived classes._  |
| virtual void | [**allocOldData**](classamrex_1_1AmrLevel.md#function-allocolddata) () <br>_Alloc space for old time data._  |
|  const BoxArray & | [**boxArray**](classamrex_1_1AmrLevel.md#function-boxarray) () noexcept const<br>_List of grids at this level._  |
| virtual void | [**checkPoint**](classamrex_1_1AmrLevel.md#function-checkpoint) (const std::string & dir, std::ostream & os, VisMF::How how=VisMF::NFiles, bool dump\_old=true) <br>_Write current state to checkpoint file._  |
| virtual void | [**checkPointPost**](classamrex_1_1AmrLevel.md#function-checkpointpost) (const std::string & dir, std::ostream & os) <br>_Do post-checkpoint work to avoid synchronizations while writing the amr hierarchy._  |
| virtual void | [**checkPointPre**](classamrex_1_1AmrLevel.md#function-checkpointpre) (const std::string & dir, std::ostream & os) <br>_Do pre-checkpoint work to avoid synchronizations while writing the amr hierarchy._  |
| virtual void | [**computeInitialDt**](classamrex_1_1AmrLevel.md#function-computeinitialdt) (int finest\_level, int sub\_cycle, Vector&lt; int &gt; & n\_cycle, const Vector&lt; IntVect &gt; & ref\_ratio, Vector&lt; Real &gt; & dt\_level, Real stop\_time) = 0<br>_Compute the initial time step. This is a pure virtual function and hence MUST be implemented by derived classes._  |
| virtual void | [**computeNewDt**](classamrex_1_1AmrLevel.md#function-computenewdt) (int finest\_level, int sub\_cycle, Vector&lt; int &gt; & n\_cycle, const Vector&lt; IntVect &gt; & ref\_ratio, Vector&lt; Real &gt; & dt\_min, Vector&lt; Real &gt; & dt\_level, Real stop\_time, int post\_regrid\_flag) = 0<br>_Compute the next time step. This is a pure virtual function and hence MUST be implemented by derived classes._  |
|  void | [**constructAreaNotToTag**](classamrex_1_1AmrLevel.md#function-constructareanottotag) () <br>_Constuct the area not to tag._  |
|  long | [**countCells**](classamrex_1_1AmrLevel.md#function-countcells) () noexcept const<br>_Returns number of cells on level._  |
| virtual std::unique\_ptr&lt; MultiFab &gt; | [**derive**](classamrex_1_1AmrLevel.md#function-derive-1-2) (const std::string & name, Real time, int ngrow) <br>_Returns a MultiFab containing the derived data for this level. The user is responsible for deleting this pointer when done with it. If ngrow&gt;0 the MultiFab is built on the appropriately grown BoxArray._  |
| virtual void | [**derive**](classamrex_1_1AmrLevel.md#function-derive-2-2) (const std::string & name, Real time, MultiFab & mf, int dcomp) <br>_This version of_ [_**derive()**_](classamrex_1_1AmrLevel.md#function-derive-1-2) _fills the dcomp'th component of mf with the derived quantity._ |
| virtual void | [**errorEst**](classamrex_1_1AmrLevel.md#function-errorest) (TagBoxArray & tb, int clearval, int tagval, Real time, int n\_error\_buf=0, int ngrow=0) = 0<br>_Error estimation for regridding. This is a pure virtual function and hence MUST be implemented by derived classes._  |
| virtual Real | [**estimateWork**](classamrex_1_1AmrLevel.md#function-estimatework) () <br>_Estimate the amount of work required to advance Just this level based on the number of cells. This estimate can be overwritten with different methods._  |
|  const IntVect & | [**fineRatio**](classamrex_1_1AmrLevel.md#function-fineratio) () noexcept const<br> |
|  const BoxArray & | [**getAreaNotToTag**](classamrex_1_1AmrLevel.md#function-getareanottotag) () noexcept<br>_Get the area not to tag._  |
|  const Box & | [**getAreaToTag**](classamrex_1_1AmrLevel.md#function-getareatotag) () noexcept<br> |
|  Vector&lt; int &gt; | [**getBCArray**](classamrex_1_1AmrLevel.md#function-getbcarray) (int State\_Type, int gridno, int scomp, int ncomp) <br>_Boundary condition access function._  |
|  const BoxArray & | [**getEdgeBoxArray**](classamrex_1_1AmrLevel.md#function-getedgeboxarray) (int dir) noexcept const<br> |
|  const BoxArray & | [**getNodalBoxArray**](classamrex_1_1AmrLevel.md#function-getnodalboxarray) () noexcept const<br> |
|  MultiFab & | [**get\_data**](classamrex_1_1AmrLevel.md#function-get-data) (int state\_indx, Real time) noexcept<br>_Get state data at specified index and time._  |
|  MultiFab & | [**get\_new\_data**](classamrex_1_1AmrLevel.md#function-get-new-data-1-2) (int state\_indx) noexcept<br>_State data at new time._  |
|  const MultiFab & | [**get\_new\_data**](classamrex_1_1AmrLevel.md#function-get-new-data-2-2) (int state\_indx) noexcept const<br>_State data at new time._  |
|  MultiFab & | [**get\_old\_data**](classamrex_1_1AmrLevel.md#function-get-old-data-1-2) (int state\_indx) noexcept<br>_State data at old time._  |
|  const MultiFab & | [**get\_old\_data**](classamrex_1_1AmrLevel.md#function-get-old-data-2-2) (int state\_indx) noexcept const<br>_State data at old time._  |
|  StateData & | [**get\_state\_data**](classamrex_1_1AmrLevel.md#function-get-state-data) (int state\_indx) noexcept<br>_State data object._  |
| virtual void | [**init**](classamrex_1_1AmrLevel.md#function-init-1-2) ([**AmrLevel**](classamrex_1_1AmrLevel.md) & old) = 0<br>_Init data on this level from another_ [_**AmrLevel**_](classamrex_1_1AmrLevel.md) _(during regrid). This is a pure virtual function and hence MUST be implemented by derived classes._ |
| virtual void | [**init**](classamrex_1_1AmrLevel.md#function-init-2-2) () = 0<br> |
| virtual void | [**initData**](classamrex_1_1AmrLevel.md#function-initdata) () = 0<br>_Init grid data at problem start-up. This is a pure virtual function and hence MUST be implemented by derived classes._  |
| virtual void | [**manual\_tags\_placement**](classamrex_1_1AmrLevel.md#function-manual-tags-placement) (TagBoxArray & tags, const Vector&lt; IntVect &gt; & bf\_lev) <br>_Called in grid\_places after other tagging routines to modify the list of tagged points. Default implementation does nothing._  |
|  int | [**nStep**](classamrex_1_1AmrLevel.md#function-nstep) () noexcept const<br>_Timestep n at this level._  |
|  int | [**numGrids**](classamrex_1_1AmrLevel.md#function-numgrids) () noexcept const<br>_Number of grids at this level._  |
|  int | [**numStates**](classamrex_1_1AmrLevel.md#function-numstates) () noexcept const<br>_Number of states at this level._  |
| virtual int | [**okToContinue**](classamrex_1_1AmrLevel.md#function-oktocontinue) () <br>_Is it ok to continue the calculation?_  |
| virtual int | [**okToRegrid**](classamrex_1_1AmrLevel.md#function-oktoregrid) () <br>_Should I regrid with this level as base level? This test is only evaluated if regrid\_int &gt; 0 and level\_count &gt;= regrid\_int as well. Defaults to true._  |
| virtual void | [**postCoarseTimeStep**](classamrex_1_1AmrLevel.md#function-postcoarsetimestep) (Real time) <br>_Contains operations to be done only after a full coarse timestep. The default implementation does nothing._  |
|  int | [**postStepRegrid**](classamrex_1_1AmrLevel.md#function-poststepregrid) () noexcept<br>_Returns whether or not we want a post-timestep regrid._  |
| virtual void | [**post\_init**](classamrex_1_1AmrLevel.md#function-post-init) (Real stop\_time) = 0<br>_Operations to be done after initialization. This is a pure virtual function and hence MUST be implemented by derived classes._  |
| virtual void | [**post\_regrid**](classamrex_1_1AmrLevel.md#function-post-regrid) (int lbase, int iteration, int new\_finest) = 0<br>_Operations to be done after regridding This is a pure virtual function and hence MUST be implemented by derived classes._  |
| virtual void | [**post\_restart**](classamrex_1_1AmrLevel.md#function-post-restart) () <br>_Operations to be done after restart._  |
| virtual void | [**post\_timestep**](classamrex_1_1AmrLevel.md#function-post-timestep) (int iteration) = 0<br>_Contains operations to be done after a timestep. This is a pure virtual function and hence MUST be implemented by derived classes._  |
| virtual void | [**removeOldData**](classamrex_1_1AmrLevel.md#function-removeolddata) () <br>_Delete old-time data._  |
|  void | [**reset**](classamrex_1_1AmrLevel.md#function-reset) () <br>_Reset data to initial time by swapping new and old time data._  |
| virtual void | [**restart**](classamrex_1_1AmrLevel.md#function-restart) ([**Amr**](classamrex_1_1Amr.md) & papa, std::istream & is, bool bReadSpecial=false) <br>_Restart from a checkpoint file._  |
|  void | [**setAreaNotToTag**](classamrex_1_1AmrLevel.md#function-setareanottotag) (BoxArray & ba) noexcept<br>_Set the area not to tag._  |
| virtual void | [**setPhysBoundaryValues**](classamrex_1_1AmrLevel.md#function-setphysboundaryvalues) (FArrayBox & dest, int state\_indx, Real time, int dest\_comp, int src\_comp, int num\_comp) <br>_Function to set physical boundary conditions._  |
| virtual void | [**setPlotVariables**](classamrex_1_1AmrLevel.md#function-setplotvariables) () <br>_Modify list of variables to be plotted._  |
|  void | [**setPostStepRegrid**](classamrex_1_1AmrLevel.md#function-setpoststepregrid) (int new\_val) noexcept<br>_Sets a new value for the post-timestep regrid trigger._  |
| virtual void | [**setSmallPlotVariables**](classamrex_1_1AmrLevel.md#function-setsmallplotvariables) () <br>_Modify list of variables to be plotted._  |
| virtual void | [**setTimeLevel**](classamrex_1_1AmrLevel.md#function-settimelevel) (Real time, Real dt\_old, Real dt\_new) <br>_Set the time levels of state data._  |
| virtual void | [**set\_preferred\_boundary\_values**](classamrex_1_1AmrLevel.md#function-set-preferred-boundary-values) (MultiFab & S, int state\_index, int scomp, int dcomp, int ncomp, Real time) const<br>_Hack to allow override of (non-fine-fine) fillpatched boundary data._  |
| virtual void | [**set\_state\_in\_checkpoint**](classamrex_1_1AmrLevel.md#function-set-state-in-checkpoint) (Vector&lt; int &gt; & state\_in\_checkpoint) <br>_Old checkpoint may have different number of states than the new source code._  |
| virtual std::string | [**thePlotFileType**](classamrex_1_1AmrLevel.md#function-theplotfiletype) () const<br>_A string written as the first item in_ [_**writePlotFile()**_](classamrex_1_1AmrLevel.md#function-writeplotfile) _at level zero. It is so we can distinguish between different types of plot files. This default "HyperCLaw-V1.1" is for VisIt software and some of our internal postprocessing routines._ |
|  [**TimeLevel**](classamrex_1_1AmrLevel.md#enum-timelevel) | [**which\_time**](classamrex_1_1AmrLevel.md#function-which-time) (int state\_indx, Real time) noexcept const<br>_Returns one the TimeLevel enums. Asserts that time is between AmrOldTime and AmrNewTime._  |
| virtual void | [**writePlotFile**](classamrex_1_1AmrLevel.md#function-writeplotfile) (const std::string & dir, std::ostream & os, VisMF::How how=VisMF::NFiles) <br>_Write plot file stuff to specified directory. This is a pure virtual function and hence MUST be implemented by derived classes._  |
| virtual void | [**writePlotFilePost**](classamrex_1_1AmrLevel.md#function-writeplotfilepost) (const std::string & dir, std::ostream & os) <br>_Do post-plotfile work to avoid synchronizations while writing the amr hierarchy._  |
| virtual void | [**writePlotFilePre**](classamrex_1_1AmrLevel.md#function-writeplotfilepre) (const std::string & dir, std::ostream & os) <br>_Do pre-plotfile work to avoid synchronizations while writing the amr hierarchy._  |
| virtual bool | [**writePlotNow**](classamrex_1_1AmrLevel.md#function-writeplotnow) () <br>_Does the_ [_**AmrLevel**_](classamrex_1_1AmrLevel.md) _want_[_**Amr**_](classamrex_1_1Amr.md) _to write a plotfile now?_ |
| virtual void | [**writeSmallPlotFile**](classamrex_1_1AmrLevel.md#function-writesmallplotfile) (const std::string & dir, std::ostream & os, VisMF::How how=VisMF::NFiles) <br>_Write small plot file stuff to specified directory._  |
| virtual bool | [**writeSmallPlotNow**](classamrex_1_1AmrLevel.md#function-writesmallplotnow) () <br>_Does the_ [_**AmrLevel**_](classamrex_1_1AmrLevel.md) _want_[_**Amr**_](classamrex_1_1Amr.md) _to write a small plotfile now?_ |
| virtual  | [**~AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel) () <br>_The destructor._  |

## Public Static Functions

| Type | Name |
| ---: | :--- |
|  void | [**FillPatch**](classamrex_1_1AmrLevel.md#function-fillpatch) ([**AmrLevel**](classamrex_1_1AmrLevel.md) & amrlevel, MultiFab & leveldata, int boxGrow, Real time, int index, int scomp, int ncomp, int dcomp=0) <br> |
|  void | [**FillPatchAdd**](classamrex_1_1AmrLevel.md#function-fillpatchadd) ([**AmrLevel**](classamrex_1_1AmrLevel.md) & amrlevel, MultiFab & leveldata, int boxGrow, Real time, int index, int scomp, int ncomp, int dcomp=0) <br> |
|  void | [**FlushFPICache**](classamrex_1_1AmrLevel.md#function-flushfpicache) () <br> |
|  DeriveList & | [**get\_derive\_lst**](classamrex_1_1AmrLevel.md#function-get-derive-lst) () noexcept<br>_Returns list of derived variables._  |
|  const DescriptorList & | [**get\_desc\_lst**](classamrex_1_1AmrLevel.md#function-get-desc-lst) () noexcept<br>_Returns list of Descriptors._  |
|  bool | [**isStateVariable**](classamrex_1_1AmrLevel.md#function-isstatevariable) (const std::string & name, int & state\_indx, int & ncomp) <br>_Is name a state variable?_  |



## Protected Attributes

| Type | Name |
| ---: | :--- |
|  IntVect | [**crse\_ratio**](classamrex_1_1AmrLevel.md#variable-crse-ratio)  <br> |
|  DistributionMapping | [**dmap**](classamrex_1_1AmrLevel.md#variable-dmap)  <br> |
|  IntVect | [**fine\_ratio**](classamrex_1_1AmrLevel.md#variable-fine-ratio)  <br> |
|  Geometry | [**geom**](classamrex_1_1AmrLevel.md#variable-geom)  <br> |
|  BoxArray | [**grids**](classamrex_1_1AmrLevel.md#variable-grids)  <br> |
|  int | [**level**](classamrex_1_1AmrLevel.md#variable-level)  <br> |
|  bool | [**levelDirectoryCreated**](classamrex_1_1AmrLevel.md#variable-leveldirectorycreated)  <br> |
|  BoxArray | [**m\_AreaNotToTag**](classamrex_1_1AmrLevel.md#variable-m-areanottotag)  <br> |
|  Box | [**m\_AreaToTag**](classamrex_1_1AmrLevel.md#variable-m-areatotag)  <br> |
|  std::unique\_ptr&lt; FabFactory&lt; FArrayBox &gt; &gt; | [**m\_factory**](classamrex_1_1AmrLevel.md#variable-m-factory)  <br> |
|  [**Amr**](classamrex_1_1Amr.md) \* | [**parent**](classamrex_1_1AmrLevel.md#variable-parent)  <br> |
|  int | [**post\_step\_regrid**](classamrex_1_1AmrLevel.md#variable-post-step-regrid)  <br> |
|  Vector&lt; StateData &gt; | [**state**](classamrex_1_1AmrLevel.md#variable-state)  <br> |

## Protected Static Attributes

| Type | Name |
| ---: | :--- |
|  DeriveList | [**derive\_lst**](classamrex_1_1AmrLevel.md#variable-derive-lst)  <br> |
|  DescriptorList | [**desc\_lst**](classamrex_1_1AmrLevel.md#variable-desc-lst)  <br> |

## Protected Functions

| Type | Name |
| ---: | :--- |
|   | [**AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel-1-3) () noexcept<br>_The constructors_  _for derived classes._ |
|   | [**AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel-2-3) ([**Amr**](classamrex_1_1Amr.md) & papa, int lev, const Geometry & level\_geom, const BoxArray & bl, const DistributionMapping & dm, Real time) <br> |
|   | [**AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel-3-3) (const [**AmrLevel**](classamrex_1_1AmrLevel.md) &) = delete<br> |
|  void | [**finishConstructor**](classamrex_1_1AmrLevel.md#function-finishconstructor) () <br>_Common code used by all constructors._  |
|  [**AmrLevel**](classamrex_1_1AmrLevel.md) & | [**operator=**](classamrex_1_1AmrLevel.md#function-operator) (const [**AmrLevel**](classamrex_1_1AmrLevel.md) &) = delete<br> |


## Public Types Documentation


### <a href="#enum-timelevel" id="enum-timelevel">enum TimeLevel </a>


```cpp
enum amrex::AmrLevel::TimeLevel {
    AmrOldTime,
    AmrHalfTime,
    AmrNewTime,
    Amr1QtrTime,
    Amr3QtrTime,
    AmrOtherTime
};
```


## Public Functions Documentation


### <a href="#function-createleveldirectory" id="function-createleveldirectory">function CreateLevelDirectory </a>


```cpp
virtual void amrex::AmrLevel::CreateLevelDirectory (
    const std::string & dir
) 
```



### <a href="#function-distributionmap" id="function-distributionmap">function DistributionMap </a>


```cpp
inline const DistributionMapping & amrex::AmrLevel::DistributionMap () noexcept const
```



### <a href="#function-domain" id="function-domain">function Domain </a>


```cpp
inline const Box & amrex::AmrLevel::Domain () noexcept const
```



### <a href="#function-factory" id="function-factory">function Factory </a>


```cpp
inline const FabFactory< FArrayBox > & amrex::AmrLevel::Factory () noexcept const
```



### <a href="#function-fillcoarsepatch" id="function-fillcoarsepatch">function FillCoarsePatch </a>


```cpp
void amrex::AmrLevel::FillCoarsePatch (
    MultiFab & dest,
    int dcomp,
    Real time,
    int state_idx,
    int scomp,
    int ncomp,
    int nghost=0
) 
```



### <a href="#function-geom" id="function-geom">function Geom </a>


```cpp
inline const Geometry & amrex::AmrLevel::Geom () noexcept const
```



### <a href="#function-level" id="function-level">function Level </a>


```cpp
inline int amrex::AmrLevel::Level () noexcept const
```



### <a href="#function-leveldirectorynames" id="function-leveldirectorynames">function LevelDirectoryNames </a>


```cpp
void amrex::AmrLevel::LevelDirectoryNames (
    const std::string & dir,
    std::string & LevelDir,
    std::string & FullPath
) 
```



### <a href="#function-setleveldirectorycreated" id="function-setleveldirectorycreated">function SetLevelDirectoryCreated </a>


```cpp
inline void amrex::AmrLevel::SetLevelDirectoryCreated (
    bool ldc
) noexcept
```



### <a href="#function-updatedistributionmaps" id="function-updatedistributionmaps">function UpdateDistributionMaps </a>


```cpp
void amrex::AmrLevel::UpdateDistributionMaps (
    DistributionMapping & dmap
) 
```



### <a href="#function-workesttype" id="function-workesttype">function WorkEstType </a>


```cpp
inline virtual int amrex::AmrLevel::WorkEstType () 
```



### <a href="#function-advance" id="function-advance">function advance </a>


```cpp
virtual Real amrex::AmrLevel::advance (
    Real time,
    Real dt,
    int iteration,
    int ncycle
) = 0
```



### <a href="#function-allocolddata" id="function-allocolddata">function allocOldData </a>


```cpp
virtual void amrex::AmrLevel::allocOldData () 
```



### <a href="#function-boxarray" id="function-boxarray">function boxArray </a>


```cpp
inline const BoxArray & amrex::AmrLevel::boxArray () noexcept const
```



### <a href="#function-checkpoint" id="function-checkpoint">function checkPoint </a>


```cpp
virtual void amrex::AmrLevel::checkPoint (
    const std::string & dir,
    std::ostream & os,
    VisMF::How how=VisMF::NFiles,
    bool dump_old=true
) 
```



### <a href="#function-checkpointpost" id="function-checkpointpost">function checkPointPost </a>


```cpp
virtual void amrex::AmrLevel::checkPointPost (
    const std::string & dir,
    std::ostream & os
) 
```



### <a href="#function-checkpointpre" id="function-checkpointpre">function checkPointPre </a>


```cpp
virtual void amrex::AmrLevel::checkPointPre (
    const std::string & dir,
    std::ostream & os
) 
```



### <a href="#function-computeinitialdt" id="function-computeinitialdt">function computeInitialDt </a>


```cpp
virtual void amrex::AmrLevel::computeInitialDt (
    int finest_level,
    int sub_cycle,
    Vector< int > & n_cycle,
    const Vector< IntVect > & ref_ratio,
    Vector< Real > & dt_level,
    Real stop_time
) = 0
```



### <a href="#function-computenewdt" id="function-computenewdt">function computeNewDt </a>


```cpp
virtual void amrex::AmrLevel::computeNewDt (
    int finest_level,
    int sub_cycle,
    Vector< int > & n_cycle,
    const Vector< IntVect > & ref_ratio,
    Vector< Real > & dt_min,
    Vector< Real > & dt_level,
    Real stop_time,
    int post_regrid_flag
) = 0
```



### <a href="#function-constructareanottotag" id="function-constructareanottotag">function constructAreaNotToTag </a>


```cpp
void amrex::AmrLevel::constructAreaNotToTag () 
```



### <a href="#function-countcells" id="function-countcells">function countCells </a>


```cpp
long amrex::AmrLevel::countCells () noexcept const
```



### <a href="#function-derive-1-2" id="function-derive-1-2">function derive [1/2]</a>


```cpp
virtual std::unique_ptr< MultiFab > amrex::AmrLevel::derive (
    const std::string & name,
    Real time,
    int ngrow
) 
```



### <a href="#function-derive-2-2" id="function-derive-2-2">function derive [2/2]</a>


```cpp
virtual void amrex::AmrLevel::derive (
    const std::string & name,
    Real time,
    MultiFab & mf,
    int dcomp
) 
```



### <a href="#function-errorest" id="function-errorest">function errorEst </a>


```cpp
virtual void amrex::AmrLevel::errorEst (
    TagBoxArray & tb,
    int clearval,
    int tagval,
    Real time,
    int n_error_buf=0,
    int ngrow=0
) = 0
```



### <a href="#function-estimatework" id="function-estimatework">function estimateWork </a>


```cpp
virtual Real amrex::AmrLevel::estimateWork () 
```



### <a href="#function-fineratio" id="function-fineratio">function fineRatio </a>


```cpp
inline const IntVect & amrex::AmrLevel::fineRatio () noexcept const
```



### <a href="#function-getareanottotag" id="function-getareanottotag">function getAreaNotToTag </a>


```cpp
const BoxArray & amrex::AmrLevel::getAreaNotToTag () noexcept
```



### <a href="#function-getareatotag" id="function-getareatotag">function getAreaToTag </a>


```cpp
const Box & amrex::AmrLevel::getAreaToTag () noexcept
```



### <a href="#function-getbcarray" id="function-getbcarray">function getBCArray </a>


```cpp
Vector< int > amrex::AmrLevel::getBCArray (
    int State_Type,
    int gridno,
    int scomp,
    int ncomp
) 
```



### <a href="#function-getedgeboxarray" id="function-getedgeboxarray">function getEdgeBoxArray </a>


```cpp
const BoxArray & amrex::AmrLevel::getEdgeBoxArray (
    int dir
) noexcept const
```



### <a href="#function-getnodalboxarray" id="function-getnodalboxarray">function getNodalBoxArray </a>


```cpp
const BoxArray & amrex::AmrLevel::getNodalBoxArray () noexcept const
```



### <a href="#function-get-data" id="function-get-data">function get\_data </a>


```cpp
MultiFab & amrex::AmrLevel::get_data (
    int state_indx,
    Real time
) noexcept
```



### <a href="#function-get-new-data-1-2" id="function-get-new-data-1-2">function get\_new\_data [1/2]</a>


```cpp
inline MultiFab & amrex::AmrLevel::get_new_data (
    int state_indx
) noexcept
```



### <a href="#function-get-new-data-2-2" id="function-get-new-data-2-2">function get\_new\_data [2/2]</a>


```cpp
inline const MultiFab & amrex::AmrLevel::get_new_data (
    int state_indx
) noexcept const
```



### <a href="#function-get-old-data-1-2" id="function-get-old-data-1-2">function get\_old\_data [1/2]</a>


```cpp
inline MultiFab & amrex::AmrLevel::get_old_data (
    int state_indx
) noexcept
```



### <a href="#function-get-old-data-2-2" id="function-get-old-data-2-2">function get\_old\_data [2/2]</a>


```cpp
inline const MultiFab & amrex::AmrLevel::get_old_data (
    int state_indx
) noexcept const
```



### <a href="#function-get-state-data" id="function-get-state-data">function get\_state\_data </a>


```cpp
inline StateData & amrex::AmrLevel::get_state_data (
    int state_indx
) noexcept
```



### <a href="#function-init-1-2" id="function-init-1-2">function init [1/2]</a>


```cpp
virtual void amrex::AmrLevel::init (
    AmrLevel & old
) = 0
```



### <a href="#function-init-2-2" id="function-init-2-2">function init [2/2]</a>


```cpp
virtual void amrex::AmrLevel::init () = 0
```


Init data on this level after regridding if old [**AmrLevel**](classamrex_1_1AmrLevel.md) did not previously exist. This is a pure virtual function and hence MUST be implemented by derived classes. 


        

### <a href="#function-initdata" id="function-initdata">function initData </a>


```cpp
virtual void amrex::AmrLevel::initData () = 0
```



### <a href="#function-manual-tags-placement" id="function-manual-tags-placement">function manual\_tags\_placement </a>


```cpp
virtual void amrex::AmrLevel::manual_tags_placement (
    TagBoxArray & tags,
    const Vector< IntVect > & bf_lev
) 
```



### <a href="#function-nstep" id="function-nstep">function nStep </a>


```cpp
inline int amrex::AmrLevel::nStep () noexcept const
```



### <a href="#function-numgrids" id="function-numgrids">function numGrids </a>


```cpp
inline int amrex::AmrLevel::numGrids () noexcept const
```



### <a href="#function-numstates" id="function-numstates">function numStates </a>


```cpp
inline int amrex::AmrLevel::numStates () noexcept const
```



### <a href="#function-oktocontinue" id="function-oktocontinue">function okToContinue </a>


```cpp
inline virtual int amrex::AmrLevel::okToContinue () 
```



### <a href="#function-oktoregrid" id="function-oktoregrid">function okToRegrid </a>


```cpp
virtual int amrex::AmrLevel::okToRegrid () 
```



### <a href="#function-postcoarsetimestep" id="function-postcoarsetimestep">function postCoarseTimeStep </a>


```cpp
virtual void amrex::AmrLevel::postCoarseTimeStep (
    Real time
) 
```



### <a href="#function-poststepregrid" id="function-poststepregrid">function postStepRegrid </a>


```cpp
inline int amrex::AmrLevel::postStepRegrid () noexcept
```



### <a href="#function-post-init" id="function-post-init">function post\_init </a>


```cpp
virtual void amrex::AmrLevel::post_init (
    Real stop_time
) = 0
```



### <a href="#function-post-regrid" id="function-post-regrid">function post\_regrid </a>


```cpp
virtual void amrex::AmrLevel::post_regrid (
    int lbase,
    int iteration,
    int new_finest
) = 0
```



### <a href="#function-post-restart" id="function-post-restart">function post\_restart </a>


```cpp
inline virtual void amrex::AmrLevel::post_restart () 
```



### <a href="#function-post-timestep" id="function-post-timestep">function post\_timestep </a>


```cpp
virtual void amrex::AmrLevel::post_timestep (
    int iteration
) = 0
```



### <a href="#function-removeolddata" id="function-removeolddata">function removeOldData </a>


```cpp
virtual void amrex::AmrLevel::removeOldData () 
```



### <a href="#function-reset" id="function-reset">function reset </a>


```cpp
void amrex::AmrLevel::reset () 
```



### <a href="#function-restart" id="function-restart">function restart </a>


```cpp
virtual void amrex::AmrLevel::restart (
    Amr & papa,
    std::istream & is,
    bool bReadSpecial=false
) 
```



### <a href="#function-setareanottotag" id="function-setareanottotag">function setAreaNotToTag </a>


```cpp
void amrex::AmrLevel::setAreaNotToTag (
    BoxArray & ba
) noexcept
```



### <a href="#function-setphysboundaryvalues" id="function-setphysboundaryvalues">function setPhysBoundaryValues </a>


```cpp
virtual void amrex::AmrLevel::setPhysBoundaryValues (
    FArrayBox & dest,
    int state_indx,
    Real time,
    int dest_comp,
    int src_comp,
    int num_comp
) 
```



### <a href="#function-setplotvariables" id="function-setplotvariables">function setPlotVariables </a>


```cpp
virtual void amrex::AmrLevel::setPlotVariables () 
```



### <a href="#function-setpoststepregrid" id="function-setpoststepregrid">function setPostStepRegrid </a>


```cpp
inline void amrex::AmrLevel::setPostStepRegrid (
    int new_val
) noexcept
```



### <a href="#function-setsmallplotvariables" id="function-setsmallplotvariables">function setSmallPlotVariables </a>


```cpp
virtual void amrex::AmrLevel::setSmallPlotVariables () 
```



### <a href="#function-settimelevel" id="function-settimelevel">function setTimeLevel </a>


```cpp
virtual void amrex::AmrLevel::setTimeLevel (
    Real time,
    Real dt_old,
    Real dt_new
) 
```



### <a href="#function-set-preferred-boundary-values" id="function-set-preferred-boundary-values">function set\_preferred\_boundary\_values </a>


```cpp
virtual void amrex::AmrLevel::set_preferred_boundary_values (
    MultiFab & S,
    int state_index,
    int scomp,
    int dcomp,
    int ncomp,
    Real time
) const
```



### <a href="#function-set-state-in-checkpoint" id="function-set-state-in-checkpoint">function set\_state\_in\_checkpoint </a>


```cpp
virtual void amrex::AmrLevel::set_state_in_checkpoint (
    Vector< int > & state_in_checkpoint
) 
```



### <a href="#function-theplotfiletype" id="function-theplotfiletype">function thePlotFileType </a>


```cpp
inline virtual std::string amrex::AmrLevel::thePlotFileType () const
```



### <a href="#function-which-time" id="function-which-time">function which\_time </a>


```cpp
TimeLevel amrex::AmrLevel::which_time (
    int state_indx,
    Real time
) noexcept const
```



### <a href="#function-writeplotfile" id="function-writeplotfile">function writePlotFile </a>


```cpp
virtual void amrex::AmrLevel::writePlotFile (
    const std::string & dir,
    std::ostream & os,
    VisMF::How how=VisMF::NFiles
) 
```



### <a href="#function-writeplotfilepost" id="function-writeplotfilepost">function writePlotFilePost </a>


```cpp
virtual void amrex::AmrLevel::writePlotFilePost (
    const std::string & dir,
    std::ostream & os
) 
```



### <a href="#function-writeplotfilepre" id="function-writeplotfilepre">function writePlotFilePre </a>


```cpp
virtual void amrex::AmrLevel::writePlotFilePre (
    const std::string & dir,
    std::ostream & os
) 
```



### <a href="#function-writeplotnow" id="function-writeplotnow">function writePlotNow </a>


```cpp
virtual bool amrex::AmrLevel::writePlotNow () 
```



### <a href="#function-writesmallplotfile" id="function-writesmallplotfile">function writeSmallPlotFile </a>


```cpp
inline virtual void amrex::AmrLevel::writeSmallPlotFile (
    const std::string & dir,
    std::ostream & os,
    VisMF::How how=VisMF::NFiles
) 
```


Unlike writePlotFile, this is NOT a pure virtual function so implementation by derived classes is optional. 


        

### <a href="#function-writesmallplotnow" id="function-writesmallplotnow">function writeSmallPlotNow </a>


```cpp
virtual bool amrex::AmrLevel::writeSmallPlotNow () 
```



### <a href="#function-amrlevel" id="function-amrlevel">function ~AmrLevel </a>


```cpp
virtual amrex::AmrLevel::~AmrLevel () 
```


## Public Static Functions Documentation


### <a href="#function-fillpatch" id="function-fillpatch">function FillPatch </a>


```cpp
static void amrex::AmrLevel::FillPatch (
    AmrLevel & amrlevel,
    MultiFab & leveldata,
    int boxGrow,
    Real time,
    int index,
    int scomp,
    int ncomp,
    int dcomp=0
) 
```



### <a href="#function-fillpatchadd" id="function-fillpatchadd">function FillPatchAdd </a>


```cpp
static void amrex::AmrLevel::FillPatchAdd (
    AmrLevel & amrlevel,
    MultiFab & leveldata,
    int boxGrow,
    Real time,
    int index,
    int scomp,
    int ncomp,
    int dcomp=0
) 
```



### <a href="#function-flushfpicache" id="function-flushfpicache">function FlushFPICache </a>


```cpp
static void amrex::AmrLevel::FlushFPICache () 
```



### <a href="#function-get-derive-lst" id="function-get-derive-lst">function get\_derive\_lst </a>


```cpp
static DeriveList & amrex::AmrLevel::get_derive_lst () noexcept
```



### <a href="#function-get-desc-lst" id="function-get-desc-lst">function get\_desc\_lst </a>


```cpp
static inline const DescriptorList & amrex::AmrLevel::get_desc_lst () noexcept
```



### <a href="#function-isstatevariable" id="function-isstatevariable">function isStateVariable </a>


```cpp
static bool amrex::AmrLevel::isStateVariable (
    const std::string & name,
    int & state_indx,
    int & ncomp
) 
```


## Protected Attributes Documentation


### <a href="#variable-crse-ratio" id="variable-crse-ratio">variable crse\_ratio </a>


```cpp
IntVect amrex::AmrLevel::crse_ratio;
```



### <a href="#variable-dmap" id="variable-dmap">variable dmap </a>


```cpp
DistributionMapping amrex::AmrLevel::dmap;
```



### <a href="#variable-fine-ratio" id="variable-fine-ratio">variable fine\_ratio </a>


```cpp
IntVect amrex::AmrLevel::fine_ratio;
```



### <a href="#variable-geom" id="variable-geom">variable geom </a>


```cpp
Geometry amrex::AmrLevel::geom;
```



### <a href="#variable-grids" id="variable-grids">variable grids </a>


```cpp
BoxArray amrex::AmrLevel::grids;
```



### <a href="#variable-level" id="variable-level">variable level </a>


```cpp
int amrex::AmrLevel::level;
```



### <a href="#variable-leveldirectorycreated" id="variable-leveldirectorycreated">variable levelDirectoryCreated </a>


```cpp
bool amrex::AmrLevel::levelDirectoryCreated;
```



### <a href="#variable-m-areanottotag" id="variable-m-areanottotag">variable m\_AreaNotToTag </a>


```cpp
BoxArray amrex::AmrLevel::m_AreaNotToTag;
```



### <a href="#variable-m-areatotag" id="variable-m-areatotag">variable m\_AreaToTag </a>


```cpp
Box amrex::AmrLevel::m_AreaToTag;
```



### <a href="#variable-m-factory" id="variable-m-factory">variable m\_factory </a>


```cpp
std::unique_ptr<FabFactory<FArrayBox> > amrex::AmrLevel::m_factory;
```



### <a href="#variable-parent" id="variable-parent">variable parent </a>


```cpp
Amr* amrex::AmrLevel::parent;
```



### <a href="#variable-post-step-regrid" id="variable-post-step-regrid">variable post\_step\_regrid </a>


```cpp
int amrex::AmrLevel::post_step_regrid;
```



### <a href="#variable-state" id="variable-state">variable state </a>


```cpp
Vector<StateData> amrex::AmrLevel::state;
```


## Protected Static Attributes Documentation


### <a href="#variable-derive-lst" id="variable-derive-lst">variable derive\_lst </a>


```cpp
DeriveList amrex::AmrLevel::derive_lst;
```



### <a href="#variable-desc-lst" id="variable-desc-lst">variable desc\_lst </a>


```cpp
DescriptorList amrex::AmrLevel::desc_lst;
```


## Protected Functions Documentation


### <a href="#function-amrlevel-1-3" id="function-amrlevel-1-3">function AmrLevel [1/3]</a>


```cpp
amrex::AmrLevel::AmrLevel () noexcept
```



### <a href="#function-amrlevel-2-3" id="function-amrlevel-2-3">function AmrLevel [2/3]</a>


```cpp
amrex::AmrLevel::AmrLevel (
    Amr & papa,
    int lev,
    const Geometry & level_geom,
    const BoxArray & bl,
    const DistributionMapping & dm,
    Real time
) 
```



### <a href="#function-amrlevel-3-3" id="function-amrlevel-3-3">function AmrLevel [3/3]</a>


```cpp
amrex::AmrLevel::AmrLevel (
    const AmrLevel &
) = delete
```



### <a href="#function-finishconstructor" id="function-finishconstructor">function finishConstructor </a>


```cpp
void amrex::AmrLevel::finishConstructor () 
```



### <a href="#function-operator" id="function-operator">function operator= </a>


```cpp
AmrLevel & amrex::AmrLevel::operator= (
    const AmrLevel &
) = delete
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/AMReX_axionyx/AMReX_AmrLevel.H`