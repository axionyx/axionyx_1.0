
# Class amrex::Amr


[**Class List**](annotated.md) **>** [**amrex**](namespaceamrex.md) **>** [**Amr**](classamrex_1_1Amr.md)



_Manage hierarchy of levels for time-dependent AMR computations._ [More...](#detailed-description)

* `#include <AMReX_Amr.H>`



Inherits the following classes: [amrex::AmrCore](classamrex_1_1AmrCore.md)












## Public Attributes

| Type | Name |
| ---: | :--- |
|  BoundaryPointList | [**intersect\_hix**](classamrex_1_1Amr.md#variable-intersect-hix)  <br> |
|  BoundaryPointList | [**intersect\_hiy**](classamrex_1_1Amr.md#variable-intersect-hiy)  <br> |
|  BoundaryPointList | [**intersect\_hiz**](classamrex_1_1Amr.md#variable-intersect-hiz)  <br> |
|  BoundaryPointList | [**intersect\_lox**](classamrex_1_1Amr.md#variable-intersect-lox)  <br> |
|  BoundaryPointList | [**intersect\_loy**](classamrex_1_1Amr.md#variable-intersect-loy)  <br> |
|  BoundaryPointList | [**intersect\_loz**](classamrex_1_1Amr.md#variable-intersect-loz)  <br> |


## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  bool | [**first\_smallplotfile**](classamrex_1_1Amr.md#variable-first-smallplotfile)  <br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Amr**](classamrex_1_1Amr.md#function-amr-1-3) () <br>_The constructor._  |
|   | [**Amr**](classamrex_1_1Amr.md#function-amr-2-3) (const RealBox \* rb, int max\_level\_in, const Vector&lt; int &gt; & n\_cell\_in, int coord) <br> |
|   | [**Amr**](classamrex_1_1Amr.md#function-amr-3-3) (const [**Amr**](classamrex_1_1Amr.md) & rhs) = delete<br> |
|  std::ostream & | [**DataLog**](classamrex_1_1Amr.md#function-datalog) (int i) <br>_The ith datalog file. Do with it what you want._  |
|  const std::string | [**DataLogName**](classamrex_1_1Amr.md#function-datalogname) (int i) noexcept const<br>_The filename of the ith datalog file._  |
|  void | [**FinalizeInit**](classamrex_1_1Amr.md#function-finalizeinit) (Real strt\_time, Real stop\_time) <br>_Second part of initialInit._  |
|  void | [**InitAmr**](classamrex_1_1Amr.md#function-initamr) () <br> |
|  void | [**InitializeInit**](classamrex_1_1Amr.md#function-initializeinit) (Real strt\_time, Real stop\_time, const BoxArray \* lev0\_grids=0, const Vector&lt; int &gt; \* pmap=0) <br>_First part of initialInit._  |
|  void | [**InstallNewDistributionMap**](classamrex_1_1Amr.md#function-installnewdistributionmap) (int lev, const DistributionMapping & newdm) <br> |
|  int | [**NumDataLogs**](classamrex_1_1Amr.md#function-numdatalogs) () noexcept<br>_How many datalogs have been opened._  |
|  bool | [**RegridOnRestart**](classamrex_1_1Amr.md#function-regridonrestart) () noexcept const<br>_Whether to regrid right after restart._  |
|  void | [**RegridOnly**](classamrex_1_1Amr.md#function-regridonly) (Real time, bool do\_io=true) <br>_Regrid only!_  |
|  bool | [**UsingPrecreateDirectories**](classamrex_1_1Amr.md#function-usingprecreatedirectories) () noexcept<br> |
|  long | [**cellCount**](classamrex_1_1Amr.md#function-cellcount-1-2) () noexcept<br>_Total number of cells._  |
|  long | [**cellCount**](classamrex_1_1Amr.md#function-cellcount-2-2) (int lev) noexcept<br>_Number of cells at given level._  |
|  int | [**checkInt**](classamrex_1_1Amr.md#function-checkint) () noexcept const<br>_Number of time steps between checkpoint files._  |
|  Real | [**checkPer**](classamrex_1_1Amr.md#function-checkper) () noexcept const<br>_Time between checkpoint files._  |
| virtual void | [**checkPoint**](classamrex_1_1Amr.md#function-checkpoint) () <br>_Write current state into a chk\* file._  |
| virtual void | [**coarseTimeStep**](classamrex_1_1Amr.md#function-coarsetimestep) (Real stop\_time) <br>_Do a complete integration cycle._  |
|  Real | [**coarseTimeStepDt**](classamrex_1_1Amr.md#function-coarsetimestepdt) (Real stop\_time) <br>_Do a complete integration cycle and return the coarse dt._  |
|  Real | [**cumTime**](classamrex_1_1Amr.md#function-cumtime) () noexcept const<br>_Physical time._  |
|  std::unique\_ptr&lt; MultiFab &gt; | [**derive**](classamrex_1_1Amr.md#function-derive) (const std::string & name, Real time, int lev, int ngrow) <br>_Retrieve derived data. User is responsible for deleting pointer._  |
|  Real | [**dtLevel**](classamrex_1_1Amr.md#function-dtlevel-1-2) (int level) noexcept const<br>_Time step at specified level._  |
|  const Vector&lt; Real &gt; & | [**dtLevel**](classamrex_1_1Amr.md#function-dtlevel-2-2) () noexcept const<br>_Array of time steps at all levels._  |
|  Real | [**dtMin**](classamrex_1_1Amr.md#function-dtmin) (int level) noexcept const<br>_Max time step (typically based on physics) at specified level._  |
|  Vector&lt; std::unique\_ptr&lt; [**AmrLevel**](classamrex_1_1AmrLevel.md) &gt; &gt; & | [**getAmrLevels**](classamrex_1_1Amr.md#function-getamrlevels) () noexcept<br>_Array of AmrLevels._  |
|  const Vector&lt; BoxArray &gt; & | [**getInitialBA**](classamrex_1_1Amr.md#function-getinitialba) () noexcept<br> |
|  BoundaryPointList & | [**getIntersectHiX**](classamrex_1_1Amr.md#function-getintersecthix) () noexcept<br> |
|  BoundaryPointList & | [**getIntersectHiY**](classamrex_1_1Amr.md#function-getintersecthiy) () noexcept<br> |
|  BoundaryPointList & | [**getIntersectHiZ**](classamrex_1_1Amr.md#function-getintersecthiz) () noexcept<br> |
|  BoundaryPointList & | [**getIntersectLoX**](classamrex_1_1Amr.md#function-getintersectlox) () noexcept<br> |
|  BoundaryPointList & | [**getIntersectLoY**](classamrex_1_1Amr.md#function-getintersectloy) () noexcept<br> |
|  BoundaryPointList & | [**getIntersectLoZ**](classamrex_1_1Amr.md#function-getintersectloz) () noexcept<br> |
|  [**AmrLevel**](classamrex_1_1AmrLevel.md) & | [**getLevel**](classamrex_1_1Amr.md#function-getlevel) (int lev) noexcept<br>[_**AmrLevel**_](classamrex_1_1AmrLevel.md) _lev._ |
| virtual void | [**init**](classamrex_1_1Amr.md#function-init) (Real strt\_time, Real stop\_time) <br>_Init data after construction. Must be called before timestepping._  |
|  int | [**levelCount**](classamrex_1_1Amr.md#function-levelcount) (int lev) noexcept const<br>_Which step are we at for the specified level?_  |
|  int | [**levelSteps**](classamrex_1_1Amr.md#function-levelsteps) (int lev) noexcept const<br>_Number of time steps at specified level._  |
|  int | [**level\_being\_advanced**](classamrex_1_1Amr.md#function-level-being-advanced) () noexcept const<br>_What is "level" in_ [_**Amr::timeStep**_](classamrex_1_1Amr.md#function-timestep) _? This is only relevant if we are still in_[_**Amr::timeStep**_](classamrex_1_1Amr.md#function-timestep) _; it is set back to -1 on leaving_[_**Amr::timeStep**_](classamrex_1_1Amr.md#function-timestep) _._ |
|  int | [**nCycle**](classamrex_1_1Amr.md#function-ncycle) (int level) noexcept const<br>_Number of subcycled time steps._  |
|  int | [**numGrids**](classamrex_1_1Amr.md#function-numgrids-1-2) () noexcept<br>_Total number of grids._  |
|  int | [**numGrids**](classamrex_1_1Amr.md#function-numgrids-2-2) (int lev) noexcept<br>_Number of grids at given level._  |
|  int | [**okToContinue**](classamrex_1_1Amr.md#function-oktocontinue) () noexcept<br>_More work to be done?_  |
|  bool | [**okToRegrid**](classamrex_1_1Amr.md#function-oktoregrid) (int level) noexcept<br>_Should we regrid this level?_  |
|  [**Amr**](classamrex_1_1Amr.md) & | [**operator=**](classamrex_1_1Amr.md#function-operator) (const [**Amr**](classamrex_1_1Amr.md) & rhs) = delete<br> |
|  int | [**plotInt**](classamrex_1_1Amr.md#function-plotint) () noexcept const<br>_Number of time steps between plot files._  |
|  Real | [**plotLogPer**](classamrex_1_1Amr.md#function-plotlogper) () noexcept const<br>_Spacing in log10(time) of logarithmically spaced plot files._  |
|  Real | [**plotPer**](classamrex_1_1Amr.md#function-plotper) () noexcept const<br>_Time between plot files._  |
|  int | [**regridInt**](classamrex_1_1Amr.md#function-regridint) (int lev) noexcept const<br>_Interval between regridding._  |
|  void | [**setBoundaryGeometry**](classamrex_1_1Amr.md#function-setboundarygeometry-1-2) (BoundaryPointList & IntersectLoX, BoundaryPointList & IntersectHiX, BoundaryPointList & IntersectLoY, BoundaryPointList & IntersectHiY) noexcept<br>_Specialized version: Define BoundaryPointLists that give the intersections of the external geometry with constant (i,k) and (j,k) These are defined at the coarsest level indexing only._  |
|  void | [**setBoundaryGeometry**](classamrex_1_1Amr.md#function-setboundarygeometry-2-2) (BoundaryPointList & IntersectLoX, BoundaryPointList & IntersectHiX, BoundaryPointList & IntersectLoY, BoundaryPointList & IntersectHiY, BoundaryPointList & IntersectLoZ, BoundaryPointList & IntersectHiZ) noexcept<br>_More general version: Define BoundaryPointLists that give the intersections of the external geometry with constant (i,k),(j,k) and (i,j). These are defined at the coarsest level indexing only._  |
|  void | [**setCumTime**](classamrex_1_1Amr.md#function-setcumtime) (Real t) noexcept<br> |
|  void | [**setDtLevel**](classamrex_1_1Amr.md#function-setdtlevel-1-2) (const Vector&lt; Real &gt; & dt\_lev) noexcept<br>_Set the timestep on each level._  |
|  void | [**setDtLevel**](classamrex_1_1Amr.md#function-setdtlevel-2-2) (Real dt, int lev) noexcept<br>_Set the timestep at one level._  |
|  void | [**setDtMin**](classamrex_1_1Amr.md#function-setdtmin) (const Vector&lt; Real &gt; & dt\_lev) noexcept<br>_Set the dtmin on each level._  |
|  void | [**setLevelCount**](classamrex_1_1Amr.md#function-setlevelcount) (int lev, int n) noexcept<br>_Which step are we at for the specified level?_  |
|  void | [**setLevelSteps**](classamrex_1_1Amr.md#function-setlevelsteps) (int lev, int n) noexcept<br>_Number of time steps at specified level._  |
|  void | [**setNCycle**](classamrex_1_1Amr.md#function-setncycle) (const Vector&lt; int &gt; & mss) noexcept<br>_Set the cycle count on each level._  |
|  void | [**setStartTime**](classamrex_1_1Amr.md#function-setstarttime) (Real t) noexcept<br> |
|  int | [**smallplotInt**](classamrex_1_1Amr.md#function-smallplotint) () noexcept const<br>_Number of time steps between small plot files._  |
|  Real | [**smallplotLogPer**](classamrex_1_1Amr.md#function-smallplotlogper) () noexcept const<br>_Spacing in log10(time) of logarithmically spaced small plot files._  |
|  Real | [**smallplotPer**](classamrex_1_1Amr.md#function-smallplotper) () noexcept const<br>_Time between plot files._  |
|  Real | [**startTime**](classamrex_1_1Amr.md#function-starttime) () noexcept const<br>_Physical time this simulation started._  |
|  int | [**stepOfLastCheckPoint**](classamrex_1_1Amr.md#function-stepoflastcheckpoint) () noexcept const<br> |
|  int | [**stepOfLastPlotFile**](classamrex_1_1Amr.md#function-stepoflastplotfile) () noexcept const<br> |
|  int | [**stepOfLastSmallPlotFile**](classamrex_1_1Amr.md#function-stepoflastsmallplotfile) () noexcept const<br> |
|  int | [**subCycle**](classamrex_1_1Amr.md#function-subcycle) () noexcept const<br>_Subcycle in time?_  |
|  const std::string & | [**subcyclingMode**](classamrex_1_1Amr.md#function-subcyclingmode) () noexcept const<br>_How are we subcycling?_  |
|  const std::string & | [**theRestartFile**](classamrex_1_1Amr.md#function-therestartfile) () noexcept const<br>_Name of the restart chkpoint file._  |
|  const std::string & | [**theRestartPlotFile**](classamrex_1_1Amr.md#function-therestartplotfile) () noexcept const<br>_Name of the restart plotfile._  |
| virtual void | [**writePlotFile**](classamrex_1_1Amr.md#function-writeplotfile) () <br>_Write the plot file to be used for visualization._  |
| virtual void | [**writeSmallPlotFile**](classamrex_1_1Amr.md#function-writesmallplotfile) () <br>_Write the small plot file to be used for visualization._  |
| virtual  | [**~Amr**](classamrex_1_1Amr.md#function-amr) () <br>_The destructor._  |

## Public Functions inherited from [amrex::AmrCore](classamrex_1_1AmrCore.md)

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
|  void | [**Finalize**](classamrex_1_1Amr.md#function-finalize) () <br> |
|  void | [**Initialize**](classamrex_1_1Amr.md#function-initialize) () <br> |
|  bool | [**Plot\_Files\_Output**](classamrex_1_1Amr.md#function-plot-files-output) () <br>_Write out plotfiles (True/False)?_  |
|  void | [**addDerivePlotVar**](classamrex_1_1Amr.md#function-addderiveplotvar) (const std::string & name) <br>_If the string is not the name of a variable in derive\_plot\_vars, add it to derive\_plot\_vars._  |
|  void | [**addDeriveSmallPlotVar**](classamrex_1_1Amr.md#function-addderivesmallplotvar) (const std::string & name) <br> |
|  void | [**addStatePlotVar**](classamrex_1_1Amr.md#function-addstateplotvar) (const std::string & name) <br>_If the string is not the name of a variable in state\_plot\_vars, add it to state\_plot\_vars._  |
|  void | [**addStateSmallPlotVar**](classamrex_1_1Amr.md#function-addstatesmallplotvar) (const std::string & name) <br> |
|  void | [**clearDerivePlotVarList**](classamrex_1_1Amr.md#function-clearderiveplotvarlist) () <br>_Clear the list of derive\_plot\_vars._  |
|  void | [**clearDeriveSmallPlotVarList**](classamrex_1_1Amr.md#function-clearderivesmallplotvarlist) () <br> |
|  void | [**clearStatePlotVarList**](classamrex_1_1Amr.md#function-clearstateplotvarlist) () <br>_Clear the list of state\_plot\_vars._  |
|  void | [**clearStateSmallPlotVarList**](classamrex_1_1Amr.md#function-clearstatesmallplotvarlist) () <br> |
|  Real | [**computeOptimalSubcycling**](classamrex_1_1Amr.md#function-computeoptimalsubcycling) (int n, int \* best, Real \* dt\_max, Real \* est\_work, int \* cycle\_max) <br>_Compute the optimal subcycling pattern. This assumes that anything less than cycle\_max[i] is a valid number of subcycles at level[i]. For example if ref\_ratio[i] = cycle\_max[i] = 4, then 1,2,3,4 are all valid values for n\_cycles[i]._  |
|  void | [**deleteDerivePlotVar**](classamrex_1_1Amr.md#function-deletederiveplotvar) (const std::string & name) <br>_Remove the string from derive\_plot\_vars._  |
|  void | [**deleteDeriveSmallPlotVar**](classamrex_1_1Amr.md#function-deletederivesmallplotvar) (const std::string & name) <br> |
|  void | [**deleteStatePlotVar**](classamrex_1_1Amr.md#function-deletestateplotvar) (const std::string & name) <br>_Remove the string from state\_plot\_vars._  |
|  const std::list&lt; std::string &gt; & | [**derivePlotVars**](classamrex_1_1Amr.md#function-deriveplotvars) () noexcept<br>_The names of derived variables to output in the plotfile. They can be set using the amr.derive\_plot\_vars variable in a ParmParse inputs file._  |
|  const std::list&lt; std::string &gt; & | [**deriveSmallPlotVars**](classamrex_1_1Amr.md#function-derivesmallplotvars) () noexcept<br> |
|  void | [**fillDerivePlotVarList**](classamrex_1_1Amr.md#function-fillderiveplotvarlist) () <br>_Fill the list of derive\_plot\_vars with all derived quantities._  |
|  void | [**fillDeriveSmallPlotVarList**](classamrex_1_1Amr.md#function-fillderivesmallplotvarlist) () <br> |
|  void | [**fillStatePlotVarList**](classamrex_1_1Amr.md#function-fillstateplotvarlist) () <br>_Fill the list of state\_plot\_vars with all of the state quantities._  |
|  void | [**fillStateSmallPlotVarList**](classamrex_1_1Amr.md#function-fillstatesmallplotvarlist) () <br> |
|  const BoxArray & | [**initialBa**](classamrex_1_1Amr.md#function-initialba) (int level) noexcept<br>_Array of BoxArrays read in to initially define grid hierarchy._  |
|  int | [**initialBaLevels**](classamrex_1_1Amr.md#function-initialbalevels) () noexcept<br>_Number of levels at which the grids are initially specified._  |
|  bool | [**isDerivePlotVar**](classamrex_1_1Amr.md#function-isderiveplotvar) (const std::string & name) noexcept<br>_Is the string the name of a variable in derive\_plot\_vars?_  |
|  bool | [**isDeriveSmallPlotVar**](classamrex_1_1Amr.md#function-isderivesmallplotvar) (const std::string & name) noexcept<br> |
|  bool | [**isStatePlotVar**](classamrex_1_1Amr.md#function-isstateplotvar) (const std::string & name) <br>_Is the string the name of a variable in state\_plot\_vars?_  |
|  bool | [**isStateSmallPlotVar**](classamrex_1_1Amr.md#function-isstatesmallplotvar) (const std::string & name) <br> |
|  const std::list&lt; std::string &gt; & | [**statePlotVars**](classamrex_1_1Amr.md#function-stateplotvars) () noexcept<br>_The names of state variables to output in the plotfile. They can be set using the amr.plot\_vars variable in a ParmParse inputs file._  |
|  const std::list&lt; std::string &gt; & | [**stateSmallPlotVars**](classamrex_1_1Amr.md#function-statesmallplotvars) () noexcept<br> |

## Public Static Functions inherited from [amrex::AmrCore](classamrex_1_1AmrCore.md)

| Type | Name |
| ---: | :--- |
|  void | [**Finalize**](classamrex_1_1AmrCore.md#function-finalize) () <br> |
|  void | [**Initialize**](classamrex_1_1AmrCore.md#function-initialize) () <br> |





## Protected Attributes

| Type | Name |
| ---: | :--- |
|  bool | [**abort\_on\_stream\_retry\_failure**](classamrex_1_1Amr.md#variable-abort-on-stream-retry-failure)  <br> |
|  Vector&lt; std::unique\_ptr&lt; [**AmrLevel**](classamrex_1_1AmrLevel.md) &gt; &gt; | [**amr\_level**](classamrex_1_1Amr.md#variable-amr-level)  <br>_Vector of levels._  |
|  bool | [**bUserStopRequest**](classamrex_1_1Amr.md#variable-buserstoprequest)  <br> |
|  std::string | [**check\_file\_root**](classamrex_1_1Amr.md#variable-check-file-root)  <br>_Root name of checkpoint file._  |
|  int | [**check\_int**](classamrex_1_1Amr.md#variable-check-int)  <br>_How often checkpoint (# time steps)._  |
|  Real | [**check\_per**](classamrex_1_1Amr.md#variable-check-per)  <br>_How often checkpoint (units of time)._  |
|  Real | [**cumtime**](classamrex_1_1Amr.md#variable-cumtime)  <br>_Physical time variable._  |
|  Vector&lt; std::unique\_ptr&lt; std::fstream &gt; &gt; | [**datalog**](classamrex_1_1Amr.md#variable-datalog)  <br> |
|  Vector&lt; std::string &gt; | [**datalogname**](classamrex_1_1Amr.md#variable-datalogname)  <br> |
|  Vector&lt; Real &gt; | [**dt\_level**](classamrex_1_1Amr.md#variable-dt-level)  <br>_Timestep at this level._  |
|  Vector&lt; Real &gt; | [**dt\_min**](classamrex_1_1Amr.md#variable-dt-min)  <br> |
|  int | [**file\_name\_digits**](classamrex_1_1Amr.md#variable-file-name-digits)  <br>_How many digits to use in the plotfile and checkpoint names._  |
|  std::ofstream | [**gridlog**](classamrex_1_1Amr.md#variable-gridlog)  <br> |
|  std::string | [**initial\_grids\_file**](classamrex_1_1Amr.md#variable-initial-grids-file)  <br>_Grids file that will bypass regridding only at initialization._  |
|  bool | [**isPeriodic**](classamrex_1_1Amr.md#variable-isperiodic)  <br>_Domain periodic?_  |
|  int | [**last\_checkpoint**](classamrex_1_1Amr.md#variable-last-checkpoint)  <br>_Step number of previous checkpoint._  |
|  int | [**last\_plotfile**](classamrex_1_1Amr.md#variable-last-plotfile)  <br>_Step number of previous plotfile._  |
|  int | [**last\_smallplotfile**](classamrex_1_1Amr.md#variable-last-smallplotfile)  <br>_Step number of previous small plotfile._  |
|  Vector&lt; int &gt; | [**level\_count**](classamrex_1_1Amr.md#variable-level-count)  <br> |
|  Vector&lt; int &gt; | [**level\_steps**](classamrex_1_1Amr.md#variable-level-steps)  <br>_Number of time steps at this level._  |
|  LevelBld \* | [**levelbld**](classamrex_1_1Amr.md#variable-levelbld)  <br> |
|  int | [**loadbalance\_level0\_int**](classamrex_1_1Amr.md#variable-loadbalance-level0-int)  <br> |
|  Real | [**loadbalance\_max\_fac**](classamrex_1_1Amr.md#variable-loadbalance-max-fac)  <br> |
|  int | [**loadbalance\_with\_workestimates**](classamrex_1_1Amr.md#variable-loadbalance-with-workestimates)  <br> |
|  int | [**message\_int**](classamrex_1_1Amr.md#variable-message-int)  <br>_How often checking messages touched by user, such as "stop\_run"._  |
|  Vector&lt; int &gt; | [**n\_cycle**](classamrex_1_1Amr.md#variable-n-cycle)  <br> |
|  std::string | [**plot\_file\_root**](classamrex_1_1Amr.md#variable-plot-file-root)  <br>_Root name of plotfile._  |
|  int | [**plot\_int**](classamrex_1_1Amr.md#variable-plot-int)  <br>_How often plotfile (# of time steps)_  |
|  Real | [**plot\_log\_per**](classamrex_1_1Amr.md#variable-plot-log-per)  <br>_How often plotfile (in units of log10(time))_  |
|  Real | [**plot\_per**](classamrex_1_1Amr.md#variable-plot-per)  <br>_How often plotfile (in units of time)_  |
|  std::string | [**probin\_file**](classamrex_1_1Amr.md#variable-probin-file)  <br> |
|  int | [**record\_grid\_info**](classamrex_1_1Amr.md#variable-record-grid-info)  <br> |
|  int | [**record\_run\_info**](classamrex_1_1Amr.md#variable-record-run-info)  <br> |
|  int | [**record\_run\_info\_terse**](classamrex_1_1Amr.md#variable-record-run-info-terse)  <br> |
|  std::string | [**regrid\_grids\_file**](classamrex_1_1Amr.md#variable-regrid-grids-file)  <br>_Grids file that will bypass regridding._  |
|  Vector&lt; int &gt; | [**regrid\_int**](classamrex_1_1Amr.md#variable-regrid-int)  <br>_Interval between regridding._  |
|  std::string | [**restart\_chkfile**](classamrex_1_1Amr.md#variable-restart-chkfile)  <br> |
|  std::string | [**restart\_pltfile**](classamrex_1_1Amr.md#variable-restart-pltfile)  <br> |
|  std::ofstream | [**runlog**](classamrex_1_1Amr.md#variable-runlog)  <br> |
|  std::ofstream | [**runlog\_terse**](classamrex_1_1Amr.md#variable-runlog-terse)  <br> |
|  std::string | [**small\_plot\_file\_root**](classamrex_1_1Amr.md#variable-small-plot-file-root)  <br>_Root name of small plotfile._  |
|  int | [**small\_plot\_int**](classamrex_1_1Amr.md#variable-small-plot-int)  <br>_How often small plotfile (# of time steps)_  |
|  Real | [**small\_plot\_log\_per**](classamrex_1_1Amr.md#variable-small-plot-log-per)  <br>_How often small plotfile (in units of log10(time))_  |
|  Real | [**small\_plot\_per**](classamrex_1_1Amr.md#variable-small-plot-per)  <br>_How often small plotfile (in units of time)_  |
|  Real | [**start\_time**](classamrex_1_1Amr.md#variable-start-time)  <br>_Physical time this simulation started._  |
|  int | [**stream\_max\_tries**](classamrex_1_1Amr.md#variable-stream-max-tries)  <br> |
|  int | [**sub\_cycle**](classamrex_1_1Amr.md#variable-sub-cycle)  <br> |
|  std::string | [**subcycling\_mode**](classamrex_1_1Amr.md#variable-subcycling-mode)  <br>_Type of subcycling to use._  |
|  int | [**which\_level\_being\_advanced**](classamrex_1_1Amr.md#variable-which-level-being-advanced)  <br>_Only &gt;=0 if we are in Amr::timeStep(level,...)_  |
|  int | [**write\_plotfile\_with\_checkpoint**](classamrex_1_1Amr.md#variable-write-plotfile-with-checkpoint)  <br>_Write out a plotfile whenever we checkpoint._  |

## Protected Attributes inherited from [amrex::AmrCore](classamrex_1_1AmrCore.md)

| Type | Name |
| ---: | :--- |
|  int | [**verbose**](classamrex_1_1AmrCore.md#variable-verbose)  <br> |

## Protected Static Attributes

| Type | Name |
| ---: | :--- |
|  std::list&lt; std::string &gt; | [**derive\_plot\_vars**](classamrex_1_1Amr.md#variable-derive-plot-vars)  <br>_Derived Vars to dump to plotfile._  |
|  std::list&lt; std::string &gt; | [**derive\_small\_plot\_vars**](classamrex_1_1Amr.md#variable-derive-small-plot-vars)  <br>_Derived Vars to dump to small plotfile._  |
|  bool | [**first\_plotfile**](classamrex_1_1Amr.md#variable-first-plotfile)  <br> |
|  Vector&lt; BoxArray &gt; | [**initial\_ba**](classamrex_1_1Amr.md#variable-initial-ba)  <br>_Array of BoxArrays read in to initially define grid hierarchy._  |
|  Vector&lt; BoxArray &gt; | [**regrid\_ba**](classamrex_1_1Amr.md#variable-regrid-ba)  <br>_Array of BoxArrays read in to externally define grid hierarchy at each regrid._  |
|  std::list&lt; std::string &gt; | [**state\_plot\_vars**](classamrex_1_1Amr.md#variable-state-plot-vars)  <br>_State Vars to dump to plotfile._  |
|  std::list&lt; std::string &gt; | [**state\_small\_plot\_vars**](classamrex_1_1Amr.md#variable-state-small-plot-vars)  <br>_State Vars to dump to small plotfile._  |


## Protected Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**ClearLevel**](classamrex_1_1Amr.md#function-clearlevel) (int lev) override<br>_Delete level data._  |
| virtual void | [**ErrorEst**](classamrex_1_1Amr.md#function-errorest) (int lev, TagBoxArray & tags, Real time, int ngrow) override<br>_Tag cells for refinement. TagBoxArray tags is built on level lev grids._  |
| virtual BoxArray | [**GetAreaNotToTag**](classamrex_1_1Amr.md#function-getareanottotag) (int lev) override<br> |
|  void | [**LoadBalanceLevel0**](classamrex_1_1Amr.md#function-loadbalancelevel0) (Real time) <br> |
| virtual void | [**MakeNewLevelFromCoarse**](classamrex_1_1Amr.md#function-makenewlevelfromcoarse) (int lev, Real time, const BoxArray & ba, const DistributionMapping & dm) override<br>_Make a new level using provided BoxArray and DistributionMapping and fill with interpolated coarse level data._  |
| virtual void | [**MakeNewLevelFromScratch**](classamrex_1_1Amr.md#function-makenewlevelfromscratch) (int lev, Real time, const BoxArray & ba, const DistributionMapping & dm) override<br> |
| virtual void | [**ManualTagsPlacement**](classamrex_1_1Amr.md#function-manualtagsplacement) (int lev, TagBoxArray & tags, const Vector&lt; IntVect &gt; & bf\_lev) override<br> |
| virtual void | [**RemakeLevel**](classamrex_1_1Amr.md#function-remakelevel) (int lev, Real time, const BoxArray & ba, const DistributionMapping & dm) override<br>_Remake an existing level using provided BoxArray and DistributionMapping and fill with existing fine and coarse data._  |
|  void | [**bldFineLevels**](classamrex_1_1Amr.md#function-bldfinelevels) (Real start\_time) <br>_Define and initialize refined levels._  |
|  void | [**checkInput**](classamrex_1_1Amr.md#function-checkinput) () <br>_Check for valid input._  |
|  void | [**defBaseLevel**](classamrex_1_1Amr.md#function-defbaselevel) (Real start\_time, const BoxArray \* lev0\_grids=0, const Vector&lt; int &gt; \* pmap=0) <br>_Define and initialize coarsest level._  |
|  void | [**grid\_places**](classamrex_1_1Amr.md#function-grid-places) (int lbase, Real time, int & new\_finest, Vector&lt; BoxArray &gt; & new\_grids) <br>_Define new grid locations (called from regrid) and put into new\_grids._  |
|  int | [**initInSitu**](classamrex_1_1Amr.md#function-initinsitu) () <br> |
|  void | [**initPltAndChk**](classamrex_1_1Amr.md#function-initpltandchk) () <br> |
|  void | [**initSubcycle**](classamrex_1_1Amr.md#function-initsubcycle) () <br> |
|  void | [**initialInit**](classamrex_1_1Amr.md#function-initialinit) (Real strt\_time, Real stop\_time, const BoxArray \* lev0\_grids=0, const Vector&lt; int &gt; \* pmap=0) <br>_Initialize grid hierarchy_  _called by_[_**Amr::init**_](classamrex_1_1Amr.md#function-init) _._ |
|  DistributionMapping | [**makeLoadBalanceDistributionMap**](classamrex_1_1Amr.md#function-makeloadbalancedistributionmap) (int lev, Real time, const BoxArray & ba) const<br> |
|  void | [**printGridInfo**](classamrex_1_1Amr.md#function-printgridinfo) (std::ostream & os, int min\_lev, int max\_lev) <br> |
|  void | [**readProbinFile**](classamrex_1_1Amr.md#function-readprobinfile) (int & init) <br>_Read the probin file._  |
| virtual void | [**regrid**](classamrex_1_1Amr.md#function-regrid) (int lbase, int iteration, Real time, bool initial=false) override<br>_Rebuild grid hierarchy finer than lbase._  |
| virtual void | [**regrid\_level\_0\_on\_restart**](classamrex_1_1Amr.md#function-regrid-level-0-on-restart) () <br>_Regrid level 0 on restart._  |
|  void | [**restart**](classamrex_1_1Amr.md#function-restart) (const std::string & filename) <br>_Restart from a checkpoint file._  |
|  void | [**setRecordDataInfo**](classamrex_1_1Amr.md#function-setrecorddatainfo) (int i, const std::string & filename) <br> |
|  void | [**setRecordGridInfo**](classamrex_1_1Amr.md#function-setrecordgridinfo) (const std::string & filename) <br> |
|  void | [**setRecordRunInfo**](classamrex_1_1Amr.md#function-setrecordruninfo) (const std::string & filename) <br> |
|  void | [**setRecordRunInfoTerse**](classamrex_1_1Amr.md#function-setrecordruninfoterse) (const std::string & filename) <br> |
| virtual void | [**timeStep**](classamrex_1_1Amr.md#function-timestep) (int level, Real time, int iteration, int niter, Real stop\_time) <br>_Do a single timestep on level L._  |
|  int | [**updateInSitu**](classamrex_1_1Amr.md#function-updateinsitu) () <br> |
|  bool | [**writePlotNow**](classamrex_1_1Amr.md#function-writeplotnow) () noexcept<br>_Whether to write a plotfile now._  |
|  bool | [**writeSmallPlotNow**](classamrex_1_1Amr.md#function-writesmallplotnow) () noexcept<br> |

## Protected Functions inherited from [amrex::AmrCore](classamrex_1_1AmrCore.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**ClearLevel**](classamrex_1_1AmrCore.md#function-clearlevel) (int lev) = 0<br>_Delete level data._  |
| virtual void | [**ErrorEst**](classamrex_1_1AmrCore.md#function-errorest) (int lev, TagBoxArray & tags, Real time, int ngrow) override = 0<br>_Tag cells for refinement. TagBoxArray tags is built on level lev grids._  |
| virtual void | [**MakeNewLevelFromCoarse**](classamrex_1_1AmrCore.md#function-makenewlevelfromcoarse) (int lev, Real time, const BoxArray & ba, const DistributionMapping & dm) = 0<br>_Make a new level using provided BoxArray and DistributionMapping and fill with interpolated coarse level data._  |
| virtual void | [**MakeNewLevelFromScratch**](classamrex_1_1AmrCore.md#function-makenewlevelfromscratch) (int lev, Real time, const BoxArray & ba, const DistributionMapping & dm) override = 0<br> |
| virtual void | [**RemakeLevel**](classamrex_1_1AmrCore.md#function-remakelevel) (int lev, Real time, const BoxArray & ba, const DistributionMapping & dm) = 0<br>_Remake an existing level using provided BoxArray and DistributionMapping and fill with existing fine and coarse data._  |

## Protected Static Functions

| Type | Name |
| ---: | :--- |
|  int | [**finalizeInSitu**](classamrex_1_1Amr.md#function-finalizeinsitu) () <br> |


## Public Attributes Documentation


### <a href="#variable-intersect-hix" id="variable-intersect-hix">variable intersect\_hix </a>


```cpp
BoundaryPointList amrex::Amr::intersect_hix;
```



### <a href="#variable-intersect-hiy" id="variable-intersect-hiy">variable intersect\_hiy </a>


```cpp
BoundaryPointList amrex::Amr::intersect_hiy;
```



### <a href="#variable-intersect-hiz" id="variable-intersect-hiz">variable intersect\_hiz </a>


```cpp
BoundaryPointList amrex::Amr::intersect_hiz;
```



### <a href="#variable-intersect-lox" id="variable-intersect-lox">variable intersect\_lox </a>


```cpp
BoundaryPointList amrex::Amr::intersect_lox;
```



### <a href="#variable-intersect-loy" id="variable-intersect-loy">variable intersect\_loy </a>


```cpp
BoundaryPointList amrex::Amr::intersect_loy;
```



### <a href="#variable-intersect-loz" id="variable-intersect-loz">variable intersect\_loz </a>


```cpp
BoundaryPointList amrex::Amr::intersect_loz;
```


## Public Static Attributes Documentation


### <a href="#variable-first-smallplotfile" id="variable-first-smallplotfile">variable first\_smallplotfile </a>


```cpp
bool amrex::Amr::first_smallplotfile;
```


## Public Functions Documentation


### <a href="#function-amr-1-3" id="function-amr-1-3">function Amr [1/3]</a>


```cpp
amrex::Amr::Amr () 
```



### <a href="#function-amr-2-3" id="function-amr-2-3">function Amr [2/3]</a>


```cpp
amrex::Amr::Amr (
    const RealBox * rb,
    int max_level_in,
    const Vector< int > & n_cell_in,
    int coord
) 
```



### <a href="#function-amr-3-3" id="function-amr-3-3">function Amr [3/3]</a>


```cpp
amrex::Amr::Amr (
    const Amr & rhs
) = delete
```



### <a href="#function-datalog" id="function-datalog">function DataLog </a>


```cpp
std::ostream & amrex::Amr::DataLog (
    int i
) 
```



### <a href="#function-datalogname" id="function-datalogname">function DataLogName </a>


```cpp
inline const std::string amrex::Amr::DataLogName (
    int i
) noexcept const
```



### <a href="#function-finalizeinit" id="function-finalizeinit">function FinalizeInit </a>


```cpp
void amrex::Amr::FinalizeInit (
    Real strt_time,
    Real stop_time
) 
```



### <a href="#function-initamr" id="function-initamr">function InitAmr </a>


```cpp
void amrex::Amr::InitAmr () 
```



### <a href="#function-initializeinit" id="function-initializeinit">function InitializeInit </a>


```cpp
void amrex::Amr::InitializeInit (
    Real strt_time,
    Real stop_time,
    const BoxArray * lev0_grids=0,
    const Vector< int > * pmap=0
) 
```



### <a href="#function-installnewdistributionmap" id="function-installnewdistributionmap">function InstallNewDistributionMap </a>


```cpp
void amrex::Amr::InstallNewDistributionMap (
    int lev,
    const DistributionMapping & newdm
) 
```



### <a href="#function-numdatalogs" id="function-numdatalogs">function NumDataLogs </a>


```cpp
int amrex::Amr::NumDataLogs () noexcept
```



### <a href="#function-regridonrestart" id="function-regridonrestart">function RegridOnRestart </a>


```cpp
bool amrex::Amr::RegridOnRestart () noexcept const
```



### <a href="#function-regridonly" id="function-regridonly">function RegridOnly </a>


```cpp
void amrex::Amr::RegridOnly (
    Real time,
    bool do_io=true
) 
```



### <a href="#function-usingprecreatedirectories" id="function-usingprecreatedirectories">function UsingPrecreateDirectories </a>


```cpp
bool amrex::Amr::UsingPrecreateDirectories () noexcept
```



### <a href="#function-cellcount-1-2" id="function-cellcount-1-2">function cellCount [1/2]</a>


```cpp
long amrex::Amr::cellCount () noexcept
```



### <a href="#function-cellcount-2-2" id="function-cellcount-2-2">function cellCount [2/2]</a>


```cpp
long amrex::Amr::cellCount (
    int lev
) noexcept
```



### <a href="#function-checkint" id="function-checkint">function checkInt </a>


```cpp
inline int amrex::Amr::checkInt () noexcept const
```



### <a href="#function-checkper" id="function-checkper">function checkPer </a>


```cpp
inline Real amrex::Amr::checkPer () noexcept const
```



### <a href="#function-checkpoint" id="function-checkpoint">function checkPoint </a>


```cpp
virtual void amrex::Amr::checkPoint () 
```



### <a href="#function-coarsetimestep" id="function-coarsetimestep">function coarseTimeStep </a>


```cpp
virtual void amrex::Amr::coarseTimeStep (
    Real stop_time
) 
```



### <a href="#function-coarsetimestepdt" id="function-coarsetimestepdt">function coarseTimeStepDt </a>


```cpp
Real amrex::Amr::coarseTimeStepDt (
    Real stop_time
) 
```



### <a href="#function-cumtime" id="function-cumtime">function cumTime </a>


```cpp
inline Real amrex::Amr::cumTime () noexcept const
```



### <a href="#function-derive" id="function-derive">function derive </a>


```cpp
std::unique_ptr< MultiFab > amrex::Amr::derive (
    const std::string & name,
    Real time,
    int lev,
    int ngrow
) 
```



### <a href="#function-dtlevel-1-2" id="function-dtlevel-1-2">function dtLevel [1/2]</a>


```cpp
inline Real amrex::Amr::dtLevel (
    int level
) noexcept const
```



### <a href="#function-dtlevel-2-2" id="function-dtlevel-2-2">function dtLevel [2/2]</a>


```cpp
inline const Vector< Real > & amrex::Amr::dtLevel () noexcept const
```



### <a href="#function-dtmin" id="function-dtmin">function dtMin </a>


```cpp
inline Real amrex::Amr::dtMin (
    int level
) noexcept const
```



### <a href="#function-getamrlevels" id="function-getamrlevels">function getAmrLevels </a>


```cpp
Vector< std::unique_ptr< AmrLevel > > & amrex::Amr::getAmrLevels () noexcept
```



### <a href="#function-getinitialba" id="function-getinitialba">function getInitialBA </a>


```cpp
const Vector< BoxArray > & amrex::Amr::getInitialBA () noexcept
```



### <a href="#function-getintersecthix" id="function-getintersecthix">function getIntersectHiX </a>


```cpp
inline BoundaryPointList & amrex::Amr::getIntersectHiX () noexcept
```



### <a href="#function-getintersecthiy" id="function-getintersecthiy">function getIntersectHiY </a>


```cpp
inline BoundaryPointList & amrex::Amr::getIntersectHiY () noexcept
```



### <a href="#function-getintersecthiz" id="function-getintersecthiz">function getIntersectHiZ </a>


```cpp
inline BoundaryPointList & amrex::Amr::getIntersectHiZ () noexcept
```



### <a href="#function-getintersectlox" id="function-getintersectlox">function getIntersectLoX </a>


```cpp
inline BoundaryPointList & amrex::Amr::getIntersectLoX () noexcept
```



### <a href="#function-getintersectloy" id="function-getintersectloy">function getIntersectLoY </a>


```cpp
inline BoundaryPointList & amrex::Amr::getIntersectLoY () noexcept
```



### <a href="#function-getintersectloz" id="function-getintersectloz">function getIntersectLoZ </a>


```cpp
inline BoundaryPointList & amrex::Amr::getIntersectLoZ () noexcept
```



### <a href="#function-getlevel" id="function-getlevel">function getLevel </a>


```cpp
inline AmrLevel & amrex::Amr::getLevel (
    int lev
) noexcept
```



### <a href="#function-init" id="function-init">function init </a>


```cpp
virtual void amrex::Amr::init (
    Real strt_time,
    Real stop_time
) 
```



### <a href="#function-levelcount" id="function-levelcount">function levelCount </a>


```cpp
inline int amrex::Amr::levelCount (
    int lev
) noexcept const
```



### <a href="#function-levelsteps" id="function-levelsteps">function levelSteps </a>


```cpp
inline int amrex::Amr::levelSteps (
    int lev
) noexcept const
```



### <a href="#function-level-being-advanced" id="function-level-being-advanced">function level\_being\_advanced </a>


```cpp
inline int amrex::Amr::level_being_advanced () noexcept const
```



### <a href="#function-ncycle" id="function-ncycle">function nCycle </a>


```cpp
inline int amrex::Amr::nCycle (
    int level
) noexcept const
```



### <a href="#function-numgrids-1-2" id="function-numgrids-1-2">function numGrids [1/2]</a>


```cpp
int amrex::Amr::numGrids () noexcept
```



### <a href="#function-numgrids-2-2" id="function-numgrids-2-2">function numGrids [2/2]</a>


```cpp
int amrex::Amr::numGrids (
    int lev
) noexcept
```



### <a href="#function-oktocontinue" id="function-oktocontinue">function okToContinue </a>


```cpp
int amrex::Amr::okToContinue () noexcept
```



### <a href="#function-oktoregrid" id="function-oktoregrid">function okToRegrid </a>


```cpp
bool amrex::Amr::okToRegrid (
    int level
) noexcept
```



### <a href="#function-operator" id="function-operator">function operator= </a>


```cpp
Amr & amrex::Amr::operator= (
    const Amr & rhs
) = delete
```



### <a href="#function-plotint" id="function-plotint">function plotInt </a>


```cpp
inline int amrex::Amr::plotInt () noexcept const
```



### <a href="#function-plotlogper" id="function-plotlogper">function plotLogPer </a>


```cpp
inline Real amrex::Amr::plotLogPer () noexcept const
```



### <a href="#function-plotper" id="function-plotper">function plotPer </a>


```cpp
inline Real amrex::Amr::plotPer () noexcept const
```



### <a href="#function-regridint" id="function-regridint">function regridInt </a>


```cpp
inline int amrex::Amr::regridInt (
    int lev
) noexcept const
```



### <a href="#function-setboundarygeometry-1-2" id="function-setboundarygeometry-1-2">function setBoundaryGeometry [1/2]</a>


```cpp
inline void amrex::Amr::setBoundaryGeometry (
    BoundaryPointList & IntersectLoX,
    BoundaryPointList & IntersectHiX,
    BoundaryPointList & IntersectLoY,
    BoundaryPointList & IntersectHiY
) noexcept
```



### <a href="#function-setboundarygeometry-2-2" id="function-setboundarygeometry-2-2">function setBoundaryGeometry [2/2]</a>


```cpp
inline void amrex::Amr::setBoundaryGeometry (
    BoundaryPointList & IntersectLoX,
    BoundaryPointList & IntersectHiX,
    BoundaryPointList & IntersectLoY,
    BoundaryPointList & IntersectHiY,
    BoundaryPointList & IntersectLoZ,
    BoundaryPointList & IntersectHiZ
) noexcept
```



### <a href="#function-setcumtime" id="function-setcumtime">function setCumTime </a>


```cpp
inline void amrex::Amr::setCumTime (
    Real t
) noexcept
```



### <a href="#function-setdtlevel-1-2" id="function-setdtlevel-1-2">function setDtLevel [1/2]</a>


```cpp
void amrex::Amr::setDtLevel (
    const Vector< Real > & dt_lev
) noexcept
```



### <a href="#function-setdtlevel-2-2" id="function-setdtlevel-2-2">function setDtLevel [2/2]</a>


```cpp
void amrex::Amr::setDtLevel (
    Real dt,
    int lev
) noexcept
```



### <a href="#function-setdtmin" id="function-setdtmin">function setDtMin </a>


```cpp
void amrex::Amr::setDtMin (
    const Vector< Real > & dt_lev
) noexcept
```



### <a href="#function-setlevelcount" id="function-setlevelcount">function setLevelCount </a>


```cpp
inline void amrex::Amr::setLevelCount (
    int lev,
    int n
) noexcept
```



### <a href="#function-setlevelsteps" id="function-setlevelsteps">function setLevelSteps </a>


```cpp
inline void amrex::Amr::setLevelSteps (
    int lev,
    int n
) noexcept
```



### <a href="#function-setncycle" id="function-setncycle">function setNCycle </a>


```cpp
void amrex::Amr::setNCycle (
    const Vector< int > & mss
) noexcept
```



### <a href="#function-setstarttime" id="function-setstarttime">function setStartTime </a>


```cpp
inline void amrex::Amr::setStartTime (
    Real t
) noexcept
```



### <a href="#function-smallplotint" id="function-smallplotint">function smallplotInt </a>


```cpp
inline int amrex::Amr::smallplotInt () noexcept const
```



### <a href="#function-smallplotlogper" id="function-smallplotlogper">function smallplotLogPer </a>


```cpp
inline Real amrex::Amr::smallplotLogPer () noexcept const
```



### <a href="#function-smallplotper" id="function-smallplotper">function smallplotPer </a>


```cpp
inline Real amrex::Amr::smallplotPer () noexcept const
```



### <a href="#function-starttime" id="function-starttime">function startTime </a>


```cpp
inline Real amrex::Amr::startTime () noexcept const
```



### <a href="#function-stepoflastcheckpoint" id="function-stepoflastcheckpoint">function stepOfLastCheckPoint </a>


```cpp
inline int amrex::Amr::stepOfLastCheckPoint () noexcept const
```



### <a href="#function-stepoflastplotfile" id="function-stepoflastplotfile">function stepOfLastPlotFile </a>


```cpp
inline int amrex::Amr::stepOfLastPlotFile () noexcept const
```



### <a href="#function-stepoflastsmallplotfile" id="function-stepoflastsmallplotfile">function stepOfLastSmallPlotFile </a>


```cpp
inline int amrex::Amr::stepOfLastSmallPlotFile () noexcept const
```



### <a href="#function-subcycle" id="function-subcycle">function subCycle </a>


```cpp
inline int amrex::Amr::subCycle () noexcept const
```



### <a href="#function-subcyclingmode" id="function-subcyclingmode">function subcyclingMode </a>


```cpp
inline const std::string & amrex::Amr::subcyclingMode () noexcept const
```



### <a href="#function-therestartfile" id="function-therestartfile">function theRestartFile </a>


```cpp
inline const std::string & amrex::Amr::theRestartFile () noexcept const
```



### <a href="#function-therestartplotfile" id="function-therestartplotfile">function theRestartPlotFile </a>


```cpp
inline const std::string & amrex::Amr::theRestartPlotFile () noexcept const
```



### <a href="#function-writeplotfile" id="function-writeplotfile">function writePlotFile </a>


```cpp
virtual void amrex::Amr::writePlotFile () 
```



### <a href="#function-writesmallplotfile" id="function-writesmallplotfile">function writeSmallPlotFile </a>


```cpp
virtual void amrex::Amr::writeSmallPlotFile () 
```



### <a href="#function-amr" id="function-amr">function ~Amr </a>


```cpp
virtual amrex::Amr::~Amr () 
```


## Public Static Functions Documentation


### <a href="#function-finalize" id="function-finalize">function Finalize </a>


```cpp
static void amrex::Amr::Finalize () 
```



### <a href="#function-initialize" id="function-initialize">function Initialize </a>


```cpp
static void amrex::Amr::Initialize () 
```



### <a href="#function-plot-files-output" id="function-plot-files-output">function Plot\_Files\_Output </a>


```cpp
static bool amrex::Amr::Plot_Files_Output () 
```



### <a href="#function-addderiveplotvar" id="function-addderiveplotvar">function addDerivePlotVar </a>


```cpp
static void amrex::Amr::addDerivePlotVar (
    const std::string & name
) 
```



### <a href="#function-addderivesmallplotvar" id="function-addderivesmallplotvar">function addDeriveSmallPlotVar </a>


```cpp
static void amrex::Amr::addDeriveSmallPlotVar (
    const std::string & name
) 
```



### <a href="#function-addstateplotvar" id="function-addstateplotvar">function addStatePlotVar </a>


```cpp
static void amrex::Amr::addStatePlotVar (
    const std::string & name
) 
```



### <a href="#function-addstatesmallplotvar" id="function-addstatesmallplotvar">function addStateSmallPlotVar </a>


```cpp
static void amrex::Amr::addStateSmallPlotVar (
    const std::string & name
) 
```



### <a href="#function-clearderiveplotvarlist" id="function-clearderiveplotvarlist">function clearDerivePlotVarList </a>


```cpp
static void amrex::Amr::clearDerivePlotVarList () 
```



### <a href="#function-clearderivesmallplotvarlist" id="function-clearderivesmallplotvarlist">function clearDeriveSmallPlotVarList </a>


```cpp
static void amrex::Amr::clearDeriveSmallPlotVarList () 
```



### <a href="#function-clearstateplotvarlist" id="function-clearstateplotvarlist">function clearStatePlotVarList </a>


```cpp
static void amrex::Amr::clearStatePlotVarList () 
```



### <a href="#function-clearstatesmallplotvarlist" id="function-clearstatesmallplotvarlist">function clearStateSmallPlotVarList </a>


```cpp
static void amrex::Amr::clearStateSmallPlotVarList () 
```



### <a href="#function-computeoptimalsubcycling" id="function-computeoptimalsubcycling">function computeOptimalSubcycling </a>


```cpp
static Real amrex::Amr::computeOptimalSubcycling (
    int n,
    int * best,
    Real * dt_max,
    Real * est_work,
    int * cycle_max
) 
```



### <a href="#function-deletederiveplotvar" id="function-deletederiveplotvar">function deleteDerivePlotVar </a>


```cpp
static void amrex::Amr::deleteDerivePlotVar (
    const std::string & name
) 
```



### <a href="#function-deletederivesmallplotvar" id="function-deletederivesmallplotvar">function deleteDeriveSmallPlotVar </a>


```cpp
static void amrex::Amr::deleteDeriveSmallPlotVar (
    const std::string & name
) 
```



### <a href="#function-deletestateplotvar" id="function-deletestateplotvar">function deleteStatePlotVar </a>


```cpp
static void amrex::Amr::deleteStatePlotVar (
    const std::string & name
) 
```



### <a href="#function-deriveplotvars" id="function-deriveplotvars">function derivePlotVars </a>


```cpp
static inline const std::list< std::string > & amrex::Amr::derivePlotVars () noexcept
```



### <a href="#function-derivesmallplotvars" id="function-derivesmallplotvars">function deriveSmallPlotVars </a>


```cpp
static inline const std::list< std::string > & amrex::Amr::deriveSmallPlotVars () noexcept
```



### <a href="#function-fillderiveplotvarlist" id="function-fillderiveplotvarlist">function fillDerivePlotVarList </a>


```cpp
static void amrex::Amr::fillDerivePlotVarList () 
```



### <a href="#function-fillderivesmallplotvarlist" id="function-fillderivesmallplotvarlist">function fillDeriveSmallPlotVarList </a>


```cpp
static void amrex::Amr::fillDeriveSmallPlotVarList () 
```



### <a href="#function-fillstateplotvarlist" id="function-fillstateplotvarlist">function fillStatePlotVarList </a>


```cpp
static void amrex::Amr::fillStatePlotVarList () 
```



### <a href="#function-fillstatesmallplotvarlist" id="function-fillstatesmallplotvarlist">function fillStateSmallPlotVarList </a>


```cpp
static void amrex::Amr::fillStateSmallPlotVarList () 
```



### <a href="#function-initialba" id="function-initialba">function initialBa </a>


```cpp
static inline const BoxArray & amrex::Amr::initialBa (
    int level
) noexcept
```



### <a href="#function-initialbalevels" id="function-initialbalevels">function initialBaLevels </a>


```cpp
static inline int amrex::Amr::initialBaLevels () noexcept
```



### <a href="#function-isderiveplotvar" id="function-isderiveplotvar">function isDerivePlotVar </a>


```cpp
static bool amrex::Amr::isDerivePlotVar (
    const std::string & name
) noexcept
```



### <a href="#function-isderivesmallplotvar" id="function-isderivesmallplotvar">function isDeriveSmallPlotVar </a>


```cpp
static bool amrex::Amr::isDeriveSmallPlotVar (
    const std::string & name
) noexcept
```



### <a href="#function-isstateplotvar" id="function-isstateplotvar">function isStatePlotVar </a>


```cpp
static bool amrex::Amr::isStatePlotVar (
    const std::string & name
) 
```



### <a href="#function-isstatesmallplotvar" id="function-isstatesmallplotvar">function isStateSmallPlotVar </a>


```cpp
static bool amrex::Amr::isStateSmallPlotVar (
    const std::string & name
) 
```



### <a href="#function-stateplotvars" id="function-stateplotvars">function statePlotVars </a>


```cpp
static inline const std::list< std::string > & amrex::Amr::statePlotVars () noexcept
```



### <a href="#function-statesmallplotvars" id="function-statesmallplotvars">function stateSmallPlotVars </a>


```cpp
static inline const std::list< std::string > & amrex::Amr::stateSmallPlotVars () noexcept
```


## Protected Attributes Documentation


### <a href="#variable-abort-on-stream-retry-failure" id="variable-abort-on-stream-retry-failure">variable abort\_on\_stream\_retry\_failure </a>


```cpp
bool amrex::Amr::abort_on_stream_retry_failure;
```



### <a href="#variable-amr-level" id="variable-amr-level">variable amr\_level </a>


```cpp
Vector<std::unique_ptr<AmrLevel> > amrex::Amr::amr_level;
```



### <a href="#variable-buserstoprequest" id="variable-buserstoprequest">variable bUserStopRequest </a>


```cpp
bool amrex::Amr::bUserStopRequest;
```



### <a href="#variable-check-file-root" id="variable-check-file-root">variable check\_file\_root </a>


```cpp
std::string amrex::Amr::check_file_root;
```



### <a href="#variable-check-int" id="variable-check-int">variable check\_int </a>


```cpp
int amrex::Amr::check_int;
```



### <a href="#variable-check-per" id="variable-check-per">variable check\_per </a>


```cpp
Real amrex::Amr::check_per;
```



### <a href="#variable-cumtime" id="variable-cumtime">variable cumtime </a>


```cpp
Real amrex::Amr::cumtime;
```



### <a href="#variable-datalog" id="variable-datalog">variable datalog </a>


```cpp
Vector<std::unique_ptr<std::fstream> > amrex::Amr::datalog;
```



### <a href="#variable-datalogname" id="variable-datalogname">variable datalogname </a>


```cpp
Vector<std::string> amrex::Amr::datalogname;
```



### <a href="#variable-dt-level" id="variable-dt-level">variable dt\_level </a>


```cpp
Vector<Real> amrex::Amr::dt_level;
```



### <a href="#variable-dt-min" id="variable-dt-min">variable dt\_min </a>


```cpp
Vector<Real> amrex::Amr::dt_min;
```



### <a href="#variable-file-name-digits" id="variable-file-name-digits">variable file\_name\_digits </a>


```cpp
int amrex::Amr::file_name_digits;
```



### <a href="#variable-gridlog" id="variable-gridlog">variable gridlog </a>


```cpp
std::ofstream amrex::Amr::gridlog;
```



### <a href="#variable-initial-grids-file" id="variable-initial-grids-file">variable initial\_grids\_file </a>


```cpp
std::string amrex::Amr::initial_grids_file;
```



### <a href="#variable-isperiodic" id="variable-isperiodic">variable isPeriodic </a>


```cpp
bool amrex::Amr::isPeriodic[AMREX_SPACEDIM];
```



### <a href="#variable-last-checkpoint" id="variable-last-checkpoint">variable last\_checkpoint </a>


```cpp
int amrex::Amr::last_checkpoint;
```



### <a href="#variable-last-plotfile" id="variable-last-plotfile">variable last\_plotfile </a>


```cpp
int amrex::Amr::last_plotfile;
```



### <a href="#variable-last-smallplotfile" id="variable-last-smallplotfile">variable last\_smallplotfile </a>


```cpp
int amrex::Amr::last_smallplotfile;
```



### <a href="#variable-level-count" id="variable-level-count">variable level\_count </a>


```cpp
Vector<int> amrex::Amr::level_count;
```



### <a href="#variable-level-steps" id="variable-level-steps">variable level\_steps </a>


```cpp
Vector<int> amrex::Amr::level_steps;
```



### <a href="#variable-levelbld" id="variable-levelbld">variable levelbld </a>


```cpp
LevelBld* amrex::Amr::levelbld;
```



### <a href="#variable-loadbalance-level0-int" id="variable-loadbalance-level0-int">variable loadbalance\_level0\_int </a>


```cpp
int amrex::Amr::loadbalance_level0_int;
```



### <a href="#variable-loadbalance-max-fac" id="variable-loadbalance-max-fac">variable loadbalance\_max\_fac </a>


```cpp
Real amrex::Amr::loadbalance_max_fac;
```



### <a href="#variable-loadbalance-with-workestimates" id="variable-loadbalance-with-workestimates">variable loadbalance\_with\_workestimates </a>


```cpp
int amrex::Amr::loadbalance_with_workestimates;
```



### <a href="#variable-message-int" id="variable-message-int">variable message\_int </a>


```cpp
int amrex::Amr::message_int;
```



### <a href="#variable-n-cycle" id="variable-n-cycle">variable n\_cycle </a>


```cpp
Vector<int> amrex::Amr::n_cycle;
```



### <a href="#variable-plot-file-root" id="variable-plot-file-root">variable plot\_file\_root </a>


```cpp
std::string amrex::Amr::plot_file_root;
```



### <a href="#variable-plot-int" id="variable-plot-int">variable plot\_int </a>


```cpp
int amrex::Amr::plot_int;
```



### <a href="#variable-plot-log-per" id="variable-plot-log-per">variable plot\_log\_per </a>


```cpp
Real amrex::Amr::plot_log_per;
```



### <a href="#variable-plot-per" id="variable-plot-per">variable plot\_per </a>


```cpp
Real amrex::Amr::plot_per;
```



### <a href="#variable-probin-file" id="variable-probin-file">variable probin\_file </a>


```cpp
std::string amrex::Amr::probin_file;
```



### <a href="#variable-record-grid-info" id="variable-record-grid-info">variable record\_grid\_info </a>


```cpp
int amrex::Amr::record_grid_info;
```



### <a href="#variable-record-run-info" id="variable-record-run-info">variable record\_run\_info </a>


```cpp
int amrex::Amr::record_run_info;
```



### <a href="#variable-record-run-info-terse" id="variable-record-run-info-terse">variable record\_run\_info\_terse </a>


```cpp
int amrex::Amr::record_run_info_terse;
```



### <a href="#variable-regrid-grids-file" id="variable-regrid-grids-file">variable regrid\_grids\_file </a>


```cpp
std::string amrex::Amr::regrid_grids_file;
```



### <a href="#variable-regrid-int" id="variable-regrid-int">variable regrid\_int </a>


```cpp
Vector<int> amrex::Amr::regrid_int;
```



### <a href="#variable-restart-chkfile" id="variable-restart-chkfile">variable restart\_chkfile </a>


```cpp
std::string amrex::Amr::restart_chkfile;
```



### <a href="#variable-restart-pltfile" id="variable-restart-pltfile">variable restart\_pltfile </a>


```cpp
std::string amrex::Amr::restart_pltfile;
```



### <a href="#variable-runlog" id="variable-runlog">variable runlog </a>


```cpp
std::ofstream amrex::Amr::runlog;
```



### <a href="#variable-runlog-terse" id="variable-runlog-terse">variable runlog\_terse </a>


```cpp
std::ofstream amrex::Amr::runlog_terse;
```



### <a href="#variable-small-plot-file-root" id="variable-small-plot-file-root">variable small\_plot\_file\_root </a>


```cpp
std::string amrex::Amr::small_plot_file_root;
```



### <a href="#variable-small-plot-int" id="variable-small-plot-int">variable small\_plot\_int </a>


```cpp
int amrex::Amr::small_plot_int;
```



### <a href="#variable-small-plot-log-per" id="variable-small-plot-log-per">variable small\_plot\_log\_per </a>


```cpp
Real amrex::Amr::small_plot_log_per;
```



### <a href="#variable-small-plot-per" id="variable-small-plot-per">variable small\_plot\_per </a>


```cpp
Real amrex::Amr::small_plot_per;
```



### <a href="#variable-start-time" id="variable-start-time">variable start\_time </a>


```cpp
Real amrex::Amr::start_time;
```



### <a href="#variable-stream-max-tries" id="variable-stream-max-tries">variable stream\_max\_tries </a>


```cpp
int amrex::Amr::stream_max_tries;
```



### <a href="#variable-sub-cycle" id="variable-sub-cycle">variable sub\_cycle </a>


```cpp
int amrex::Amr::sub_cycle;
```



### <a href="#variable-subcycling-mode" id="variable-subcycling-mode">variable subcycling\_mode </a>


```cpp
std::string amrex::Amr::subcycling_mode;
```



### <a href="#variable-which-level-being-advanced" id="variable-which-level-being-advanced">variable which\_level\_being\_advanced </a>


```cpp
int amrex::Amr::which_level_being_advanced;
```



### <a href="#variable-write-plotfile-with-checkpoint" id="variable-write-plotfile-with-checkpoint">variable write\_plotfile\_with\_checkpoint </a>


```cpp
int amrex::Amr::write_plotfile_with_checkpoint;
```


## Protected Static Attributes Documentation


### <a href="#variable-derive-plot-vars" id="variable-derive-plot-vars">variable derive\_plot\_vars </a>


```cpp
std::list< std::string > amrex::Amr::derive_plot_vars;
```



### <a href="#variable-derive-small-plot-vars" id="variable-derive-small-plot-vars">variable derive\_small\_plot\_vars </a>


```cpp
std::list< std::string > amrex::Amr::derive_small_plot_vars;
```



### <a href="#variable-first-plotfile" id="variable-first-plotfile">variable first\_plotfile </a>


```cpp
bool amrex::Amr::first_plotfile;
```



### <a href="#variable-initial-ba" id="variable-initial-ba">variable initial\_ba </a>


```cpp
Vector< BoxArray > amrex::Amr::initial_ba;
```



### <a href="#variable-regrid-ba" id="variable-regrid-ba">variable regrid\_ba </a>


```cpp
Vector< BoxArray > amrex::Amr::regrid_ba;
```



### <a href="#variable-state-plot-vars" id="variable-state-plot-vars">variable state\_plot\_vars </a>


```cpp
std::list< std::string > amrex::Amr::state_plot_vars;
```



### <a href="#variable-state-small-plot-vars" id="variable-state-small-plot-vars">variable state\_small\_plot\_vars </a>


```cpp
std::list< std::string > amrex::Amr::state_small_plot_vars;
```


## Protected Functions Documentation


### <a href="#function-clearlevel" id="function-clearlevel">function ClearLevel </a>


```cpp
inline virtual void amrex::Amr::ClearLevel (
    int lev
) override
```


Implements [*amrex::AmrCore::ClearLevel*](classamrex_1_1AmrCore.md#function-clearlevel)


### <a href="#function-errorest" id="function-errorest">function ErrorEst </a>


```cpp
virtual void amrex::Amr::ErrorEst (
    int lev,
    TagBoxArray & tags,
    Real time,
    int ngrow
) override
```


Implements [*amrex::AmrCore::ErrorEst*](classamrex_1_1AmrCore.md#function-errorest)


### <a href="#function-getareanottotag" id="function-getareanottotag">function GetAreaNotToTag </a>


```cpp
virtual BoxArray amrex::Amr::GetAreaNotToTag (
    int lev
) override
```



### <a href="#function-loadbalancelevel0" id="function-loadbalancelevel0">function LoadBalanceLevel0 </a>


```cpp
void amrex::Amr::LoadBalanceLevel0 (
    Real time
) 
```



### <a href="#function-makenewlevelfromcoarse" id="function-makenewlevelfromcoarse">function MakeNewLevelFromCoarse </a>


```cpp
inline virtual void amrex::Amr::MakeNewLevelFromCoarse (
    int lev,
    Real time,
    const BoxArray & ba,
    const DistributionMapping & dm
) override
```


Implements [*amrex::AmrCore::MakeNewLevelFromCoarse*](classamrex_1_1AmrCore.md#function-makenewlevelfromcoarse)


### <a href="#function-makenewlevelfromscratch" id="function-makenewlevelfromscratch">function MakeNewLevelFromScratch </a>


```cpp
inline virtual void amrex::Amr::MakeNewLevelFromScratch (
    int lev,
    Real time,
    const BoxArray & ba,
    const DistributionMapping & dm
) override
```


Make a new level from scratch using provided BoxArray and DistributionMapping. Only used during initialization. 


        
Implements [*amrex::AmrCore::MakeNewLevelFromScratch*](classamrex_1_1AmrCore.md#function-makenewlevelfromscratch)


### <a href="#function-manualtagsplacement" id="function-manualtagsplacement">function ManualTagsPlacement </a>


```cpp
virtual void amrex::Amr::ManualTagsPlacement (
    int lev,
    TagBoxArray & tags,
    const Vector< IntVect > & bf_lev
) override
```



### <a href="#function-remakelevel" id="function-remakelevel">function RemakeLevel </a>


```cpp
inline virtual void amrex::Amr::RemakeLevel (
    int lev,
    Real time,
    const BoxArray & ba,
    const DistributionMapping & dm
) override
```


Implements [*amrex::AmrCore::RemakeLevel*](classamrex_1_1AmrCore.md#function-remakelevel)


### <a href="#function-bldfinelevels" id="function-bldfinelevels">function bldFineLevels </a>


```cpp
void amrex::Amr::bldFineLevels (
    Real start_time
) 
```



### <a href="#function-checkinput" id="function-checkinput">function checkInput </a>


```cpp
void amrex::Amr::checkInput () 
```



### <a href="#function-defbaselevel" id="function-defbaselevel">function defBaseLevel </a>


```cpp
void amrex::Amr::defBaseLevel (
    Real start_time,
    const BoxArray * lev0_grids=0,
    const Vector< int > * pmap=0
) 
```



### <a href="#function-grid-places" id="function-grid-places">function grid\_places </a>


```cpp
void amrex::Amr::grid_places (
    int lbase,
    Real time,
    int & new_finest,
    Vector< BoxArray > & new_grids
) 
```



### <a href="#function-initinsitu" id="function-initinsitu">function initInSitu </a>


```cpp
int amrex::Amr::initInSitu () 
```



### <a href="#function-initpltandchk" id="function-initpltandchk">function initPltAndChk </a>


```cpp
void amrex::Amr::initPltAndChk () 
```



### <a href="#function-initsubcycle" id="function-initsubcycle">function initSubcycle </a>


```cpp
void amrex::Amr::initSubcycle () 
```



### <a href="#function-initialinit" id="function-initialinit">function initialInit </a>


```cpp
void amrex::Amr::initialInit (
    Real strt_time,
    Real stop_time,
    const BoxArray * lev0_grids=0,
    const Vector< int > * pmap=0
) 
```



### <a href="#function-makeloadbalancedistributionmap" id="function-makeloadbalancedistributionmap">function makeLoadBalanceDistributionMap </a>


```cpp
DistributionMapping amrex::Amr::makeLoadBalanceDistributionMap (
    int lev,
    Real time,
    const BoxArray & ba
) const
```



### <a href="#function-printgridinfo" id="function-printgridinfo">function printGridInfo </a>


```cpp
void amrex::Amr::printGridInfo (
    std::ostream & os,
    int min_lev,
    int max_lev
) 
```



### <a href="#function-readprobinfile" id="function-readprobinfile">function readProbinFile </a>


```cpp
void amrex::Amr::readProbinFile (
    int & init
) 
```



### <a href="#function-regrid" id="function-regrid">function regrid </a>


```cpp
virtual void amrex::Amr::regrid (
    int lbase,
    int iteration,
    Real time,
    bool initial=false
) override
```


Implements [*amrex::AmrCore::regrid*](classamrex_1_1AmrCore.md#function-regrid)


### <a href="#function-regrid-level-0-on-restart" id="function-regrid-level-0-on-restart">function regrid\_level\_0\_on\_restart </a>


```cpp
virtual void amrex::Amr::regrid_level_0_on_restart () 
```



### <a href="#function-restart" id="function-restart">function restart </a>


```cpp
void amrex::Amr::restart (
    const std::string & filename
) 
```



### <a href="#function-setrecorddatainfo" id="function-setrecorddatainfo">function setRecordDataInfo </a>


```cpp
void amrex::Amr::setRecordDataInfo (
    int i,
    const std::string & filename
) 
```



### <a href="#function-setrecordgridinfo" id="function-setrecordgridinfo">function setRecordGridInfo </a>


```cpp
void amrex::Amr::setRecordGridInfo (
    const std::string & filename
) 
```



### <a href="#function-setrecordruninfo" id="function-setrecordruninfo">function setRecordRunInfo </a>


```cpp
void amrex::Amr::setRecordRunInfo (
    const std::string & filename
) 
```



### <a href="#function-setrecordruninfoterse" id="function-setrecordruninfoterse">function setRecordRunInfoTerse </a>


```cpp
void amrex::Amr::setRecordRunInfoTerse (
    const std::string & filename
) 
```



### <a href="#function-timestep" id="function-timestep">function timeStep </a>


```cpp
virtual void amrex::Amr::timeStep (
    int level,
    Real time,
    int iteration,
    int niter,
    Real stop_time
) 
```



### <a href="#function-updateinsitu" id="function-updateinsitu">function updateInSitu </a>


```cpp
int amrex::Amr::updateInSitu () 
```



### <a href="#function-writeplotnow" id="function-writeplotnow">function writePlotNow </a>


```cpp
bool amrex::Amr::writePlotNow () noexcept
```



### <a href="#function-writesmallplotnow" id="function-writesmallplotnow">function writeSmallPlotNow </a>


```cpp
bool amrex::Amr::writeSmallPlotNow () noexcept
```


## Protected Static Functions Documentation


### <a href="#function-finalizeinsitu" id="function-finalizeinsitu">function finalizeInSitu </a>


```cpp
static int amrex::Amr::finalizeInSitu () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/AMReX_axionyx/AMReX_Amr.H`