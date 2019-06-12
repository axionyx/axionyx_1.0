
# Class Nyx


[**Class List**](annotated.md) **>** [**Nyx**](classNyx.md)



_AmrLevel-derived class for hyperbolic conservation equations for stellar media._ 

* `#include <Nyx.H>`



Inherits the following classes: [amrex::AmrLevel](classamrex_1_1AmrLevel.md),  [amrex::AmrLevel](classamrex_1_1AmrLevel.md)









## Public Types inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

| Type | Name |
| ---: | :--- |
| enum  | [**TimeLevel**](classamrex_1_1AmrLevel.md#enum-timelevel)  <br>_What time are we at?_  |

## Public Types inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

| Type | Name |
| ---: | :--- |
| enum  | [**TimeLevel**](classamrex_1_1AmrLevel.md#enum-timelevel)  <br>_What time are we at?_  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  amrex::MultiFab \* | [**fine\_mask**](classNyx.md#variable-fine-mask)  <br> |



## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  int | [**Density**](classNyx.md#variable-density)   = = -1<br> |
|  int | [**Eden**](classNyx.md#variable-eden)   = = -1<br> |
|  int | [**Eint**](classNyx.md#variable-eint)   = = -1<br> |
|  int | [**FirstAdv**](classNyx.md#variable-firstadv)   = = -1<br> |
|  int | [**FirstAux**](classNyx.md#variable-firstaux)   = = -1<br> |
|  int | [**FirstSpec**](classNyx.md#variable-firstspec)   = = -1<br> |
|  int | [**NUM\_STATE**](classNyx.md#variable-num-state)   = = -1<br> |
|  int | [**Ne\_comp**](classNyx.md#variable-ne-comp)   = = -1<br> |
|  int | [**NumAdv**](classNyx.md#variable-numadv)   = = 0<br> |
|  int | [**NumAux**](classNyx.md#variable-numaux)   = = 0<br> |
|  int | [**NumSpec**](classNyx.md#variable-numspec)   = = 0<br> |
|  int | [**Temp\_comp**](classNyx.md#variable-temp-comp)   = = -1<br> |
|  int | [**Xmom**](classNyx.md#variable-xmom)   = = -1<br> |
|  int | [**Ymom**](classNyx.md#variable-ymom)   = = -1<br> |
|  int | [**Zhi\_comp**](classNyx.md#variable-zhi-comp)   = = -1<br> |
|  int | [**Zmom**](classNyx.md#variable-zmom)   = = -1<br> |
|  amrex::Real | [**absolute\_max\_change\_a**](classNyx.md#variable-absolute-max-change-a)   = = -1.0<br>_Absolute change in a allowed in one timestep for fixed delta\_a._  |
|  amrex::Real | [**comoving\_OmB**](classNyx.md#variable-comoving-omb)  <br>_comoving parameters_  |
|  amrex::Real | [**comoving\_OmM**](classNyx.md#variable-comoving-omm)  <br> |
|  amrex::Real | [**comoving\_h**](classNyx.md#variable-comoving-h)  <br> |
|  amrex::Real | [**dt\_binpow**](classNyx.md#variable-dt-binpow)   = = -1.0<br>_Positive number means use powers of 2 binning for relative dt._  |
|  amrex::Real | [**final\_a**](classNyx.md#variable-final-a)   = = -1.0<br>_Final a_  _used as stopping criterion if positive._ |
|  amrex::Real | [**final\_time**](classNyx.md#variable-final-time)   = = -1.0<br>_End time in code units._  |
|  amrex::Real | [**final\_z**](classNyx.md#variable-final-z)   = = -1.0<br>_Final z_  _used as stopping criterion if positive._ |
|  int | [**init\_with\_sph\_particles**](classNyx.md#variable-init-with-sph-particles)   = = 0<br> |
|  amrex::Real | [**initial\_time**](classNyx.md#variable-initial-time)   = = -1.0<br>_Initial time in code units._  |
|  amrex::Real | [**initial\_z**](classNyx.md#variable-initial-z)   = = -1.0<br>_Initial redshift._  |
|  amrex::Real | [**new\_a**](classNyx.md#variable-new-a)   = = -1.0<br> |
|  amrex::Real | [**new\_a\_time**](classNyx.md#variable-new-a-time)   = = -1.0<br> |
|  amrex::Real | [**old\_a**](classNyx.md#variable-old-a)   = = -1.0<br>_"a" at old\_a\_time and new\_a\_time_  |
|  amrex::Real | [**old\_a\_time**](classNyx.md#variable-old-a-time)   = = -1.0<br>_Old and new times at which "old\_a" and "new\_a" are defined._  |
|  amrex::Real | [**particle\_cfl**](classNyx.md#variable-particle-cfl)   = = 0.5<br>_Default cfl of particles in Particle class._  |
|  std::string | [**particle\_plotfile\_format**](classNyx.md#variable-particle-plotfile-format)   = = "NATIVE"<br> |
|  int | [**particle\_verbose**](classNyx.md#variable-particle-verbose)   = = 1<br>_Default verbosity of Particle class._  |
|  int | [**print\_fortran\_warnings**](classNyx.md#variable-print-fortran-warnings)   = = true<br>_If true then print the warnings from the Fortran routines._  |
|  amrex::Real | [**relative\_max\_change\_a**](classNyx.md#variable-relative-max-change-a)   = =  0.01<br>_Relative change in a allowed in one timestep._  |
|  int | [**strict\_subcycling**](classNyx.md#variable-strict-subcycling)   = = 0<br> |
|  int | [**write\_coarsened\_particles**](classNyx.md#variable-write-coarsened-particles)   = = 0<br> |
|  int | [**write\_parameters\_in\_plotfile**](classNyx.md#variable-write-parameters-in-plotfile)   = = true<br>_Write all parameters into specified directory._  |
|  int | [**write\_particle\_density\_at\_init**](classNyx.md#variable-write-particle-density-at-init)   = = 0<br> |



## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**CreateLevelDirectory**](classNyx.md#function-createleveldirectory-1-2) (const std::string & dir) <br>_Create the Level\_ and other directories in checkpoint and plot files._  |
| virtual void | [**CreateLevelDirectory**](classNyx.md#function-createleveldirectory-2-2) (const std::string & dir) <br>_Create the Level\_ and other directories in checkpoint and plot files._  |
|  void | [**LevelDirectoryNames**](classNyx.md#function-leveldirectorynames-1-2) (const std::string & dir, const std::string & secondDir, std::string & LevelDir, std::string & FullPath) <br>_Get the level directory names._  |
|  void | [**LevelDirectoryNames**](classNyx.md#function-leveldirectorynames-1-2) (const std::string & dir, const std::string & secondDir, std::string & LevelDir, std::string & FullPath) <br>_Get the level directory names._  |
|  void | [**Lya\_statistics**](classNyx.md#function-lya-statistics-1-2) () <br> |
|  void | [**Lya\_statistics**](classNyx.md#function-lya-statistics-1-2) () <br> |
|   | [**Nyx**](classNyx.md#function-nyx-1-4) () <br>_Default constructor. Builds invalid object._  |
|   | [**Nyx**](classNyx.md#function-nyx-2-4) ([**amrex::Amr**](classamrex_1_1Amr.md) & papa, int lev, const amrex::Geometry & level\_geom, const amrex::BoxArray & bl, const amrex::DistributionMapping & dm, amrex::Real time) <br>_The basic constructor._  |
|   | [**Nyx**](classNyx.md#function-nyx-1-4) () <br>_Default constructor. Builds invalid object._  |
|   | [**Nyx**](classNyx.md#function-nyx-2-4) ([**amrex::Amr**](classamrex_1_1Amr.md) & papa, int lev, const amrex::Geometry & level\_geom, const amrex::BoxArray & bl, const amrex::DistributionMapping & dm, amrex::Real time) <br>_The basic constructor._  |
|  void | [**ReadPlotFile**](classNyx.md#function-readplotfile-1-2) (bool first, const std::string & plot\_file\_name, bool & rhoe\_infile) <br>_Initialize grid data from a plotfile at problem start-up._  |
|  void | [**ReadPlotFile**](classNyx.md#function-readplotfile-1-2) (bool first, const std::string & plot\_file\_name, bool & rhoe\_infile) <br>_Initialize grid data from a plotfile at problem start-up._  |
| virtual amrex::Real | [**advance**](classNyx.md#function-advance-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br>_Advance grids at this level in time._  |
| virtual amrex::Real | [**advance**](classNyx.md#function-advance-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br>_Advance grids at this level in time._  |
|  amrex::Real | [**advance\_hydro**](classNyx.md#function-advance-hydro-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br> |
|  amrex::Real | [**advance\_hydro**](classNyx.md#function-advance-hydro-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br> |
|  amrex::Real | [**advance\_hydro\_plus\_particles**](classNyx.md#function-advance-hydro-plus-particles-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br> |
|  amrex::Real | [**advance\_hydro\_plus\_particles**](classNyx.md#function-advance-hydro-plus-particles-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br> |
|  amrex::Real | [**advance\_no\_hydro**](classNyx.md#function-advance-no-hydro-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br> |
|  amrex::Real | [**advance\_no\_hydro**](classNyx.md#function-advance-no-hydro-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br> |
|  amrex::Real | [**advance\_particles\_only**](classNyx.md#function-advance-particles-only-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br> |
|  amrex::Real | [**advance\_particles\_only**](classNyx.md#function-advance-particles-only-1-2) (amrex::Real time, amrex::Real dt, int iteration, int ncycle) <br> |
|  void | [**analysis\_z\_est\_time\_step**](classNyx.md#function-analysis-z-est-time-step-1-2) (amrex::Real & est\_dt, bool & dt\_changed) <br>_Time step control based on "z" not passing one of the specified analysis\_z\_values._  |
|  void | [**analysis\_z\_est\_time\_step**](classNyx.md#function-analysis-z-est-time-step-1-2) (amrex::Real & est\_dt, bool & dt\_changed) <br>_Time step control based on "z" not passing one of the specified analysis\_z\_values._  |
|  amrex::MultiFab \* | [**build\_fine\_mask**](classNyx.md#function-build-fine-mask-1-2) () <br> |
|  amrex::MultiFab \* | [**build\_fine\_mask**](classNyx.md#function-build-fine-mask-2-2) () <br> |
| virtual void | [**checkPoint**](classNyx.md#function-checkpoint-1-2) (const std::string & dir, std::ostream & os, amrex::VisMF::How how, bool dump\_old) <br>_Call_ [_**amrex::AmrLevel::checkPoint**_](classamrex_1_1AmrLevel.md#function-checkpoint) _and then add radiation info._ |
| virtual void | [**checkPoint**](classNyx.md#function-checkpoint-1-2) (const std::string & dir, std::ostream & os, amrex::VisMF::How how, bool dump\_old) <br>_Call_ [_**amrex::AmrLevel::checkPoint**_](classamrex_1_1AmrLevel.md#function-checkpoint) _and then add radiation info._ |
| virtual void | [**checkPointPost**](classNyx.md#function-checkpointpost-1-2) (const std::string & dir, std::ostream & os) <br>_Do post-checkpoint work to avoid synchronizations while writing the amr hierarchy._  |
| virtual void | [**checkPointPost**](classNyx.md#function-checkpointpost-2-2) (const std::string & dir, std::ostream & os) <br>_Do post-checkpoint work to avoid synchronizations while writing the amr hierarchy._  |
| virtual void | [**checkPointPre**](classNyx.md#function-checkpointpre-1-2) (const std::string & dir, std::ostream & os) <br>_Do pre-checkpoint work to avoid synchronizations while writing the amr hierarchy._  |
| virtual void | [**checkPointPre**](classNyx.md#function-checkpointpre-2-2) (const std::string & dir, std::ostream & os) <br>_Do pre-checkpoint work to avoid synchronizations while writing the amr hierarchy._  |
|  void | [**comoving\_a\_post\_restart**](classNyx.md#function-comoving-a-post-restart-1-2) (const std::string & restart\_file) <br>_How to initialize "a" at restart (from checkpoint or plotfile)_  |
|  void | [**comoving\_a\_post\_restart**](classNyx.md#function-comoving-a-post-restart-1-2) (const std::string & restart\_file) <br>_How to initialize "a" at restart (from checkpoint or plotfile)_  |
|  void | [**comoving\_est\_time\_step**](classNyx.md#function-comoving-est-time-step-1-2) (amrex::Real & cur\_time, amrex::Real & est\_dt) <br>_Time step control based on "a" not growing too fast._  |
|  void | [**comoving\_est\_time\_step**](classNyx.md#function-comoving-est-time-step-1-2) (amrex::Real & cur\_time, amrex::Real & est\_dt) <br>_Time step control based on "a" not growing too fast._  |
| virtual void | [**computeInitialDt**](classNyx.md#function-computeinitialdt-1-2) (int finest\_level, int sub\_cycle, amrex::Vector&lt; int &gt; & n\_cycle, const amrex::Vector&lt; amrex::IntVect &gt; & ref\_ratio, amrex::Vector&lt; amrex::Real &gt; & dt\_level, amrex::Real stop\_time) <br>_Compute initial_  _dt_ _._ |
| virtual void | [**computeInitialDt**](classNyx.md#function-computeinitialdt-1-2) (int finest\_level, int sub\_cycle, amrex::Vector&lt; int &gt; & n\_cycle, const amrex::Vector&lt; amrex::IntVect &gt; & ref\_ratio, amrex::Vector&lt; amrex::Real &gt; & dt\_level, amrex::Real stop\_time) <br>_Compute initial_  _dt_ _._ |
| virtual void | [**computeNewDt**](classNyx.md#function-computenewdt-1-2) (int finest\_level, int sub\_cycle, amrex::Vector&lt; int &gt; & n\_cycle, const amrex::Vector&lt; amrex::IntVect &gt; & ref\_ratio, amrex::Vector&lt; amrex::Real &gt; & dt\_min, amrex::Vector&lt; amrex::Real &gt; & dt\_level, amrex::Real stop\_time, int post\_regrid\_flag) <br>_Compute new_  _dt_ _._ |
| virtual void | [**computeNewDt**](classNyx.md#function-computenewdt-1-2) (int finest\_level, int sub\_cycle, amrex::Vector&lt; int &gt; & n\_cycle, const amrex::Vector&lt; amrex::IntVect &gt; & ref\_ratio, amrex::Vector&lt; amrex::Real &gt; & dt\_min, amrex::Vector&lt; amrex::Real &gt; & dt\_level, amrex::Real stop\_time, int post\_regrid\_flag) <br>_Compute new_  _dt_ _._ |
|  void | [**compute\_gas\_fractions**](classNyx.md#function-compute-gas-fractions-1-2) (amrex::Real T\_cut, amrex::Real rho\_cut, amrex::Real & whim\_mass\_frac, amrex::Real & whim\_vol\_frac, amrex::Real & hh\_mass\_frac, amrex::Real & hh\_vol\_frac, amrex::Real & igm\_mass\_frac, amrex::Real & igm\_vol\_frac) <br> |
|  void | [**compute\_gas\_fractions**](classNyx.md#function-compute-gas-fractions-1-2) (amrex::Real T\_cut, amrex::Real rho\_cut, amrex::Real & whim\_mass\_frac, amrex::Real & whim\_vol\_frac, amrex::Real & hh\_mass\_frac, amrex::Real & hh\_vol\_frac, amrex::Real & igm\_mass\_frac, amrex::Real & igm\_vol\_frac) <br> |
|  void | [**compute\_hydro\_sources**](classNyx.md#function-compute-hydro-sources-1-2) (amrex::Real time, amrex::Real dt, amrex::Real a\_old, amrex::Real a\_new, amrex::MultiFab & S\_border, amrex::MultiFab & D\_border, amrex::MultiFab & ext\_src\_old, amrex::MultiFab & hydro\_src, amrex::MultiFab & grav, amrex::MultiFab & divu\_cc, bool init\_flux\_register, bool add\_to\_flux\_register) <br> |
|  void | [**compute\_hydro\_sources**](classNyx.md#function-compute-hydro-sources-1-2) (amrex::Real time, amrex::Real dt, amrex::Real a\_old, amrex::Real a\_new, amrex::MultiFab & S\_border, amrex::MultiFab & D\_border, amrex::MultiFab & ext\_src\_old, amrex::MultiFab & hydro\_src, amrex::MultiFab & grav, amrex::MultiFab & divu\_cc, bool init\_flux\_register, bool add\_to\_flux\_register) <br> |
|  void | [**compute\_new\_temp**](classNyx.md#function-compute-new-temp-1-2) (amrex::MultiFab & S\_new, amrex::MultiFab & D\_new) <br>_Note: this no longer includes the call to reset\_internal\_energy._  |
|  void | [**compute\_new\_temp**](classNyx.md#function-compute-new-temp-1-2) (amrex::MultiFab & S\_new, amrex::MultiFab & D\_new) <br>_Note: this no longer includes the call to reset\_internal\_energy._  |
|  void | [**compute\_rho\_temp**](classNyx.md#function-compute-rho-temp-1-2) (amrex::Real & rho\_T\_avg, amrex::Real & T\_avg, amrex::Real & Tinv\_avg, amrex::Real & T\_meanrho) <br> |
|  void | [**compute\_rho\_temp**](classNyx.md#function-compute-rho-temp-1-2) (amrex::Real & rho\_T\_avg, amrex::Real & T\_avg, amrex::Real & Tinv\_avg, amrex::Real & T\_meanrho) <br> |
|  void | [**conserved\_to\_primitive**](classNyx.md#function-conserved-to-primitive-1-2) (amrex::MultiFab & state) <br> |
|  void | [**conserved\_to\_primitive**](classNyx.md#function-conserved-to-primitive-1-2) (amrex::MultiFab & state) <br> |
|  std::unique\_ptr&lt; amrex::MultiFab &gt; | [**derive**](classNyx.md#function-derive-1-4) (const std::string & name, amrex::Real time, int ngrow) <br> |
|  void | [**derive**](classNyx.md#function-derive-2-4) (const std::string & name, amrex::Real time, amrex::MultiFab & mf, int dcomp) <br> |
|  std::unique\_ptr&lt; amrex::MultiFab &gt; | [**derive**](classNyx.md#function-derive-1-4) (const std::string & name, amrex::Real time, int ngrow) <br> |
|  void | [**derive**](classNyx.md#function-derive-2-4) (const std::string & name, amrex::Real time, amrex::MultiFab & mf, int dcomp) <br> |
|  bool | [**doAnalysisNow**](classNyx.md#function-doanalysisnow-1-2) () <br>_Tell_ [_**Nyx**_](classNyx.md) _to do analysis now._ |
|  bool | [**doAnalysisNow**](classNyx.md#function-doanalysisnow-1-2) () <br>_Tell_ [_**Nyx**_](classNyx.md) _to do analysis now._ |
|  void | [**do\_energy\_diagnostics**](classNyx.md#function-do-energy-diagnostics-1-2) () <br>_Print information about energy budget._  |
|  void | [**do\_energy\_diagnostics**](classNyx.md#function-do-energy-diagnostics-1-2) () <br>_Print information about energy budget._  |
| virtual void | [**errorEst**](classNyx.md#function-errorest-1-2) (amrex::TagBoxArray & tb, int clearval, int tagval, amrex::Real time, int n\_error\_buf=0, int ngrow=0) <br>_Error estimation for regridding._  |
| virtual void | [**errorEst**](classNyx.md#function-errorest-1-2) (amrex::TagBoxArray & tb, int clearval, int tagval, amrex::Real time, int n\_error\_buf=0, int ngrow=0) <br>_Error estimation for regridding._  |
|  amrex::Real | [**est\_time\_step**](classNyx.md#function-est-time-step-1-2) (amrex::Real dt\_old) <br>_Estimate time step._  |
|  amrex::Real | [**est\_time\_step**](classNyx.md#function-est-time-step-1-2) (amrex::Real dt\_old) <br>_Estimate time step._  |
|  amrex::Real | [**get\_comoving\_a**](classNyx.md#function-get-comoving-a-1-2) (amrex::Real time) <br>_Get the comoving coordinate "a"._  |
|  amrex::Real | [**get\_comoving\_a**](classNyx.md#function-get-comoving-a-1-2) (amrex::Real time) <br>_Get the comoving coordinate "a"._  |
|  void | [**get\_new\_source**](classNyx.md#function-get-new-source-1-2) (amrex::Real old\_time, amrex::Real new\_time, amrex::Real dt, amrex::MultiFab & Rhs) <br> |
|  void | [**get\_new\_source**](classNyx.md#function-get-new-source-1-2) (amrex::Real old\_time, amrex::Real new\_time, amrex::Real dt, amrex::MultiFab & Rhs) <br> |
|  void | [**get\_old\_source**](classNyx.md#function-get-old-source-1-2) (amrex::Real old\_time, amrex::Real dt, amrex::MultiFab & Rhs) <br> |
|  void | [**get\_old\_source**](classNyx.md#function-get-old-source-1-2) (amrex::Real old\_time, amrex::Real dt, amrex::MultiFab & Rhs) <br> |
|  void | [**halo\_accrete**](classNyx.md#function-halo-accrete-1-2) (amrex::Real dt) <br> |
|  void | [**halo\_accrete**](classNyx.md#function-halo-accrete-1-2) (amrex::Real dt) <br> |
|  void | [**halo\_find**](classNyx.md#function-halo-find-1-2) (amrex::Real dt) <br> |
|  void | [**halo\_find**](classNyx.md#function-halo-find-1-2) (amrex::Real dt) <br> |
|  void | [**halo\_merge**](classNyx.md#function-halo-merge-1-2) () <br> |
|  void | [**halo\_merge**](classNyx.md#function-halo-merge-1-2) () <br> |
| virtual void | [**init**](classNyx.md#function-init-1-4) ([**amrex::AmrLevel**](classamrex_1_1AmrLevel.md) & old) <br>_Initialize data on this level from another_ [_**Nyx**_](classNyx.md) _(during regrid)._ |
| virtual void | [**init**](classNyx.md#function-init-2-4) () <br>_Initialize data on this level after regridding if old level did not._  |
| virtual void | [**init**](classNyx.md#function-init-3-4) ([**amrex::AmrLevel**](classamrex_1_1AmrLevel.md) & old) <br>_Initialize data on this level from another_ [_**Nyx**_](classNyx.md) _(during regrid)._ |
| virtual void | [**init**](classNyx.md#function-init-4-4) () <br>_Initialize data on this level after regridding if old level did not._  |
| virtual void | [**initData**](classNyx.md#function-initdata-1-2) () <br>_Initialize grid data at problem start-up._  |
| virtual void | [**initData**](classNyx.md#function-initdata-2-2) () <br>_Initialize grid data at problem start-up._  |
|  void | [**init\_from\_plotfile**](classNyx.md#function-init-from-plotfile-1-2) () <br>_Initialize grid data from a plotfile at problem start-up._  |
|  void | [**init\_from\_plotfile**](classNyx.md#function-init-from-plotfile-1-2) () <br>_Initialize grid data from a plotfile at problem start-up._  |
| virtual void | [**init\_particles**](classNyx.md#function-init-particles-1-2) () <br>_Initialize particle locations and velocities (and strengths if relevant)_  |
| virtual void | [**init\_particles**](classNyx.md#function-init-particles-2-2) () <br>_Initialize particle locations and velocities (and strengths if relevant)_  |
|  void | [**init\_zhi**](classNyx.md#function-init-zhi-1-2) () <br>_Initialize the zhi component of the EOS from a binary input file._  |
|  void | [**init\_zhi**](classNyx.md#function-init-zhi-1-2) () <br>_Initialize the zhi component of the EOS from a binary input file._  |
|  void | [**initcosmo**](classNyx.md#function-initcosmo-1-2) () <br>_Initialize from MUSIC._  |
|  void | [**initcosmo**](classNyx.md#function-initcosmo-1-2) () <br>_Initialize from MUSIC._  |
|  amrex::Real | [**initial\_time\_step**](classNyx.md#function-initial-time-step-1-2) () <br>_Compute initial time step._  |
|  amrex::Real | [**initial\_time\_step**](classNyx.md#function-initial-time-step-2-2) () <br>_Compute initial time step._  |
|  void | [**integrate\_comoving\_a**](classNyx.md#function-integrate-comoving-a-1-2) (amrex::Real time, amrex::Real dt) <br> |
|  void | [**integrate\_comoving\_a**](classNyx.md#function-integrate-comoving-a-1-2) (amrex::Real time, amrex::Real dt) <br> |
|  int | [**integrate\_state\_box**](classNyx.md#function-integrate-state-box-1-2) (amrex::MultiFab & state, amrex::MultiFab & diag\_eos, const amrex::Real & a, const amrex::Real & delta\_time) <br> |
|  int | [**integrate\_state\_box**](classNyx.md#function-integrate-state-box-1-2) (amrex::MultiFab & state, amrex::MultiFab & diag\_eos, const amrex::Real & a, const amrex::Real & delta\_time) <br> |
|  int | [**integrate\_state\_grownbox**](classNyx.md#function-integrate-state-grownbox-1-2) (amrex::MultiFab & state, amrex::MultiFab & diag\_eos, const amrex::Real & a, const amrex::Real & delta\_time) <br> |
|  int | [**integrate\_state\_grownbox**](classNyx.md#function-integrate-state-grownbox-1-2) (amrex::MultiFab & state, amrex::MultiFab & diag\_eos, const amrex::Real & a, const amrex::Real & delta\_time) <br> |
|  int | [**integrate\_state\_grownvec**](classNyx.md#function-integrate-state-grownvec-1-2) (amrex::MultiFab & state, amrex::MultiFab & diag\_eos, const amrex::Real & a, const amrex::Real & delta\_time) <br> |
|  int | [**integrate\_state\_grownvec**](classNyx.md#function-integrate-state-grownvec-1-2) (amrex::MultiFab & state, amrex::MultiFab & diag\_eos, const amrex::Real & a, const amrex::Real & delta\_time) <br> |
|  int | [**integrate\_state\_vec**](classNyx.md#function-integrate-state-vec-1-2) (amrex::MultiFab & state, amrex::MultiFab & diag\_eos, const amrex::Real & a, const amrex::Real & delta\_time) <br> |
|  int | [**integrate\_state\_vec**](classNyx.md#function-integrate-state-vec-1-2) (amrex::MultiFab & state, amrex::MultiFab & diag\_eos, const amrex::Real & a, const amrex::Real & delta\_time) <br> |
| virtual void | [**manual\_tags\_placement**](classNyx.md#function-manual-tags-placement-1-2) (amrex::TagBoxArray & tags, const amrex::Vector&lt; amrex::IntVect &gt; & bf\_lev) override<br> |
| virtual void | [**manual\_tags\_placement**](classNyx.md#function-manual-tags-placement-1-2) (amrex::TagBoxArray & tags, const amrex::Vector&lt; amrex::IntVect &gt; & bf\_lev) override<br> |
|  void | [**moveKickDriftExact**](classNyx.md#function-movekickdriftexact-1-2) (amrex::Real dt) <br> |
|  void | [**moveKickDriftExact**](classNyx.md#function-movekickdriftexact-1-2) (amrex::Real dt) <br> |
|  void | [**moveKickExact**](classNyx.md#function-movekickexact-1-2) (amrex::Real dt) <br> |
|  void | [**moveKickExact**](classNyx.md#function-movekickexact-1-2) (amrex::Real dt) <br> |
| virtual int | [**okToContinue**](classNyx.md#function-oktocontinue-1-2) () <br>_Proceed with next timestep?_  |
| virtual int | [**okToContinue**](classNyx.md#function-oktocontinue-2-2) () <br>_Proceed with next timestep?_  |
|  void | [**particle\_check\_point**](classNyx.md#function-particle-check-point-1-2) (const std::string & dir) <br>_Write particles in checkpoint directories._  |
|  void | [**particle\_check\_point**](classNyx.md#function-particle-check-point-1-2) (const std::string & dir) <br>_Write particles in checkpoint directories._  |
|  std::unique\_ptr&lt; amrex::MultiFab &gt; | [**particle\_derive**](classNyx.md#function-particle-derive-1-2) (const std::string & name, amrex::Real time, int ngrow) <br>_Derived quantities associated with particles._  |
|  std::unique\_ptr&lt; amrex::MultiFab &gt; | [**particle\_derive**](classNyx.md#function-particle-derive-1-2) (const std::string & name, amrex::Real time, int ngrow) <br>_Derived quantities associated with particles._  |
|  void | [**particle\_est\_time\_step**](classNyx.md#function-particle-est-time-step-1-2) (amrex::Real & est\_dt) <br>_Time step control based on particles._  |
|  void | [**particle\_est\_time\_step**](classNyx.md#function-particle-est-time-step-1-2) (amrex::Real & est\_dt) <br>_Time step control based on particles._  |
|  void | [**particle\_move\_random**](classNyx.md#function-particle-move-random-1-2) () <br>_Move randomly._  |
|  void | [**particle\_move\_random**](classNyx.md#function-particle-move-random-1-2) () <br>_Move randomly._  |
|  void | [**particle\_plot\_file**](classNyx.md#function-particle-plot-file-1-2) (const std::string & dir) <br>_Write particles in plotfile directories._  |
|  void | [**particle\_plot\_file**](classNyx.md#function-particle-plot-file-1-2) (const std::string & dir) <br>_Write particles in plotfile directories._  |
|  void | [**particle\_post\_restart**](classNyx.md#function-particle-post-restart-1-2) (const std::string & restart\_file, bool is\_checkpoint=true) <br>_How to initialize at restart._  |
|  void | [**particle\_post\_restart**](classNyx.md#function-particle-post-restart-1-2) (const std::string & restart\_file, bool is\_checkpoint=true) <br>_How to initialize at restart._  |
|  void | [**particle\_redistribute**](classNyx.md#function-particle-redistribute-1-2) (int lbase=0, bool init=false) <br>_Redistribute._  |
|  void | [**particle\_redistribute**](classNyx.md#function-particle-redistribute-2-2) (int lbase=0, int iteration=0, bool init=false) <br>_Redistribute._  |
|  void | [**plot\_z\_est\_time\_step**](classNyx.md#function-plot-z-est-time-step-1-2) (amrex::Real & est\_dt, bool & dt\_changed) <br>_Time step control based on "z" not passing one of the specified plot\_z\_values._  |
|  void | [**plot\_z\_est\_time\_step**](classNyx.md#function-plot-z-est-time-step-1-2) (amrex::Real & est\_dt, bool & dt\_changed) <br>_Time step control based on "z" not passing one of the specified plot\_z\_values._  |
| virtual void | [**postCoarseTimeStep**](classNyx.md#function-postcoarsetimestep-1-2) (amrex::Real cumtime) <br>_Contains operations to be done only after a full coarse timestep._  |
| virtual void | [**postCoarseTimeStep**](classNyx.md#function-postcoarsetimestep-1-2) (amrex::Real cumtime) <br>_Contains operations to be done only after a full coarse timestep._  |
| virtual void | [**post\_init**](classNyx.md#function-post-init-1-2) (amrex::Real stop\_time) <br>_Do work after_ `init()` _._ |
| virtual void | [**post\_init**](classNyx.md#function-post-init-1-2) (amrex::Real stop\_time) <br>_Do work after_ `init()` _._ |
| virtual void | [**post\_regrid**](classNyx.md#function-post-regrid-1-2) (int lbase, int new\_finest) <br>_Do work after_ `regrid()` _._ |
| virtual void | [**post\_regrid**](classNyx.md#function-post-regrid-2-2) (int lbase, int iteration, int new\_finest) <br>_Do work after_ `regrid()` _._ |
| virtual void | [**post\_restart**](classNyx.md#function-post-restart-1-2) () <br>_Do work after a_ `restart()` _._ |
| virtual void | [**post\_restart**](classNyx.md#function-post-restart-2-2) () <br>_Do work after a_ `restart()` _._ |
| virtual void | [**post\_timestep**](classNyx.md#function-post-timestep-1-2) (int iteration) <br>_Do work after timestep()._  |
| virtual void | [**post\_timestep**](classNyx.md#function-post-timestep-2-2) (int iteration) <br>_Do work after timestep()._  |
|  void | [**primitive\_to\_conserved**](classNyx.md#function-primitive-to-conserved-1-2) (amrex::MultiFab & state) <br> |
|  void | [**primitive\_to\_conserved**](classNyx.md#function-primitive-to-conserved-1-2) (amrex::MultiFab & state) <br> |
|  void | [**remove\_ghost\_particles**](classNyx.md#function-remove-ghost-particles-1-2) () <br>_Remove ghost particles (for this level) if necessary._  |
|  void | [**remove\_ghost\_particles**](classNyx.md#function-remove-ghost-particles-1-2) () <br>_Remove ghost particles (for this level) if necessary._  |
|  void | [**remove\_virtual\_particles**](classNyx.md#function-remove-virtual-particles-1-2) () <br>_Remove virtual particles if necessary._  |
|  void | [**remove\_virtual\_particles**](classNyx.md#function-remove-virtual-particles-1-2) () <br>_Remove virtual particles if necessary._  |
|  void | [**reset\_internal\_energy**](classNyx.md#function-reset-internal-energy-1-2) (amrex::MultiFab & State, amrex::MultiFab & DiagEOS, amrex::MultiFab & reset\_e\_src) <br>_Synchronize (rho e) and (rho E) so they are consistent with each other._  |
|  void | [**reset\_internal\_energy**](classNyx.md#function-reset-internal-energy-1-2) (amrex::MultiFab & State, amrex::MultiFab & DiagEOS, amrex::MultiFab & reset\_e\_src) <br>_Synchronize (rho e) and (rho E) so they are consistent with each other._  |
| virtual void | [**restart**](classNyx.md#function-restart-1-2) ([**amrex::Amr**](classamrex_1_1Amr.md) & papa, istream & is, bool b\_read\_special=false) <br>_Restart from a checkpoint file._  |
| virtual void | [**restart**](classNyx.md#function-restart-2-2) ([**amrex::Amr**](classamrex_1_1Amr.md) & papa, istream & is, bool b\_read\_special=false) <br>_Restart from a checkpoint file._  |
| virtual void | [**setPlotVariables**](classNyx.md#function-setplotvariables-1-2) () <br>_Modify list of variables to be plotted._  |
| virtual void | [**setPlotVariables**](classNyx.md#function-setplotvariables-2-2) () <br>_Modify list of variables to be plotted._  |
| virtual void | [**setTimeLevel**](classNyx.md#function-settimelevel-1-2) (amrex::Real time, amrex::Real dt\_old, amrex::Real dt\_new) <br>_Set time levels of state data._  |
| virtual void | [**setTimeLevel**](classNyx.md#function-settimelevel-1-2) (amrex::Real time, amrex::Real dt\_old, amrex::Real dt\_new) <br>_Set time levels of state data._  |
|  void | [**setup\_ghost\_particles**](classNyx.md#function-setup-ghost-particles-1-2) (int ngrow) <br>_Setup ghost particles (for finer levels) if necessary._  |
|  void | [**setup\_ghost\_particles**](classNyx.md#function-setup-ghost-particles-1-2) (int ngrow) <br>_Setup ghost particles (for finer levels) if necessary._  |
|  void | [**setup\_virtual\_particles**](classNyx.md#function-setup-virtual-particles-1-2) () <br>_Setup virtual particles if necessary._  |
|  void | [**setup\_virtual\_particles**](classNyx.md#function-setup-virtual-particles-1-2) () <br>_Setup virtual particles if necessary._  |
|  void | [**strang\_first\_step**](classNyx.md#function-strang-first-step-1-2) (amrex::Real time, amrex::Real dt, amrex::MultiFab & state, amrex::MultiFab & dstate) <br> |
|  void | [**strang\_first\_step**](classNyx.md#function-strang-first-step-1-2) (amrex::Real time, amrex::Real dt, amrex::MultiFab & state, amrex::MultiFab & dstate) <br> |
|  void | [**strang\_hydro**](classNyx.md#function-strang-hydro-1-2) (amrex::Real time, amrex::Real dt, amrex::Real a\_old, amrex::Real a\_new) <br> |
|  void | [**strang\_hydro**](classNyx.md#function-strang-hydro-1-2) (amrex::Real time, amrex::Real dt, amrex::Real a\_old, amrex::Real a\_new) <br> |
|  void | [**strang\_second\_step**](classNyx.md#function-strang-second-step-1-2) (amrex::Real time, amrex::Real dt, amrex::MultiFab & state, amrex::MultiFab & dstate) <br> |
|  void | [**strang\_second\_step**](classNyx.md#function-strang-second-step-1-2) (amrex::Real time, amrex::Real dt, amrex::MultiFab & state, amrex::MultiFab & dstate) <br> |
| virtual std::string | [**thePlotFileType**](classNyx.md#function-theplotfiletype-1-2) () const<br> |
| virtual std::string | [**thePlotFileType**](classNyx.md#function-theplotfiletype-2-2) () const<br> |
|  void | [**time\_center\_source\_terms**](classNyx.md#function-time-center-source-terms-1-2) (amrex::MultiFab & S\_new, amrex::MultiFab & ext\_src\_old, amrex::MultiFab & ext\_src\_new, amrex::Real dt) <br> |
|  void | [**time\_center\_source\_terms**](classNyx.md#function-time-center-source-terms-1-2) (amrex::MultiFab & S\_new, amrex::MultiFab & ext\_src\_old, amrex::MultiFab & ext\_src\_new, amrex::Real dt) <br> |
|  void | [**update\_state\_with\_sources**](classNyx.md#function-update-state-with-sources-1-2) (amrex::MultiFab & S\_old, amrex::MultiFab & S\_new, amrex::MultiFab & ext\_src\_old, amrex::MultiFab & hydro\_src, amrex::MultiFab & grav, amrex::MultiFab & divu\_cc, amrex::Real dt, amrex::Real a\_old, amrex::Real a\_new) <br> |
|  void | [**update\_state\_with\_sources**](classNyx.md#function-update-state-with-sources-1-2) (amrex::MultiFab & S\_old, amrex::MultiFab & S\_new, amrex::MultiFab & ext\_src\_old, amrex::MultiFab & hydro\_src, amrex::MultiFab & grav, amrex::MultiFab & divu\_cc, amrex::Real dt, amrex::Real a\_old, amrex::Real a\_new) <br> |
|  amrex::Real | [**vol\_weight\_squared\_sum**](classNyx.md#function-vol-weight-squared-sum-1-2) (const std::string & name, amrex::Real time) <br> |
|  amrex::Real | [**vol\_weight\_squared\_sum**](classNyx.md#function-vol-weight-squared-sum-1-2) (const std::string & name, amrex::Real time) <br> |
|  amrex::Real | [**vol\_weight\_squared\_sum\_level**](classNyx.md#function-vol-weight-squared-sum-level-1-2) (const std::string & name, amrex::Real time) <br> |
|  amrex::Real | [**vol\_weight\_squared\_sum\_level**](classNyx.md#function-vol-weight-squared-sum-level-1-2) (const std::string & name, amrex::Real time) <br> |
|  amrex::Real | [**vol\_weight\_sum**](classNyx.md#function-vol-weight-sum-1-4) (const std::string & name, amrex::Real time, bool masked) <br> |
|  amrex::Real | [**vol\_weight\_sum**](classNyx.md#function-vol-weight-sum-2-4) (amrex::MultiFab & mf, bool masked) <br> |
|  amrex::Real | [**vol\_weight\_sum**](classNyx.md#function-vol-weight-sum-1-4) (const std::string & name, amrex::Real time, bool masked) <br> |
|  amrex::Real | [**vol\_weight\_sum**](classNyx.md#function-vol-weight-sum-2-4) (amrex::MultiFab & mf, bool masked) <br> |
|  void | [**writeJobInfo**](classNyx.md#function-writejobinfo-1-2) (const std::string & dir) <br> |
|  void | [**writeJobInfo**](classNyx.md#function-writejobinfo-1-2) (const std::string & dir) <br> |
|  void | [**writeMultiFabAsPlotFile**](classNyx.md#function-writemultifabasplotfile-1-2) (const std::string & pltfile, const amrex::MultiFab & mf, std::string componentName) <br>_Write amrex::MultiFab as plot file._  |
|  void | [**writeMultiFabAsPlotFile**](classNyx.md#function-writemultifabasplotfile-1-2) (const std::string & pltfile, const amrex::MultiFab & mf, std::string componentName) <br>_Write amrex::MultiFab as plot file._  |
| virtual void | [**writePlotFile**](classNyx.md#function-writeplotfile-1-2) (const std::string & dir, ostream & os, amrex::VisMF::How how) <br>_Write a plotfile to specified directory._  |
| virtual void | [**writePlotFile**](classNyx.md#function-writeplotfile-1-2) (const std::string & dir, ostream & os, amrex::VisMF::How how) <br>_Write a plotfile to specified directory._  |
| virtual void | [**writePlotFilePost**](classNyx.md#function-writeplotfilepost-1-2) (const std::string & dir, ostream & os) <br> |
| virtual void | [**writePlotFilePost**](classNyx.md#function-writeplotfilepost-2-2) (const std::string & dir, ostream & os) <br> |
| virtual void | [**writePlotFilePre**](classNyx.md#function-writeplotfilepre-1-2) (const std::string & dir, ostream & os) <br> |
| virtual void | [**writePlotFilePre**](classNyx.md#function-writeplotfilepre-2-2) (const std::string & dir, ostream & os) <br> |
| virtual bool | [**writePlotNow**](classNyx.md#function-writeplotnow-1-2) () <br>_Tell_ [_**amrex::Amr**_](classamrex_1_1Amr.md) _to write a plotfile now._ |
| virtual bool | [**writePlotNow**](classNyx.md#function-writeplotnow-1-2) () <br>_Tell_ [_**amrex::Amr**_](classamrex_1_1Amr.md) _to write a plotfile now._ |
| virtual void | [**write\_parameter\_file**](classNyx.md#function-write-parameter-file-1-2) (const std::string & dir) <br> |
| virtual void | [**write\_parameter\_file**](classNyx.md#function-write-parameter-file-2-2) (const std::string & dir) <br> |
| virtual  | [**~Nyx**](classNyx.md#function-nyx-1-2) () <br>_The destructor._  |
| virtual  | [**~Nyx**](classNyx.md#function-nyx-2-2) () <br>_The destructor._  |

## Public Functions inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

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

## Public Functions inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

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
|  int | [**Do\_Hydro**](classNyx.md#function-do-hydro-1-2) () <br> |
|  int | [**Do\_Hydro**](classNyx.md#function-do-hydro-2-2) () <br> |
|  void | [**InitDeriveList**](classNyx.md#function-initderivelist-1-2) () <br> |
|  void | [**InitDeriveList**](classNyx.md#function-initderivelist-2-2) () <br> |
|  void | [**InitErrorList**](classNyx.md#function-initerrorlist-1-2) () <br> |
|  void | [**InitErrorList**](classNyx.md#function-initerrorlist-2-2) () <br> |
|  void | [**alloc\_simd\_vec**](classNyx.md#function-alloc-simd-vec-1-2) () <br> |
|  void | [**alloc\_simd\_vec**](classNyx.md#function-alloc-simd-vec-1-2) () <br> |
|  void | [**dealloc\_simd\_vec**](classNyx.md#function-dealloc-simd-vec-1-2) () <br> |
|  void | [**dealloc\_simd\_vec**](classNyx.md#function-dealloc-simd-vec-1-2) () <br> |
|  void | [**error\_setup**](classNyx.md#function-error-setup-1-2) () <br>_Define tagging functions._  |
|  void | [**error\_setup**](classNyx.md#function-error-setup-2-2) () <br>_Define tagging functions._  |
|  void | [**hydro\_setup**](classNyx.md#function-hydro-setup-1-2) () <br> |
|  void | [**hydro\_setup**](classNyx.md#function-hydro-setup-2-2) () <br> |
|  void | [**no\_hydro\_setup**](classNyx.md#function-no-hydro-setup-1-2) () <br> |
|  void | [**no\_hydro\_setup**](classNyx.md#function-no-hydro-setup-1-2) () <br> |
|  int | [**num\_grow**](classNyx.md#function-num-grow-1-2) () <br> |
|  int | [**num\_grow**](classNyx.md#function-num-grow-2-2) () <br> |
|  void | [**read\_comoving\_params**](classNyx.md#function-read-comoving-params-1-2) () <br>_Read inputs related to comoving coordinates._  |
|  void | [**read\_comoving\_params**](classNyx.md#function-read-comoving-params-2-2) () <br>_Read inputs related to comoving coordinates._  |
|  void | [**read\_init\_params**](classNyx.md#function-read-init-params-1-2) () <br>_Read initialization-related inputs._  |
|  void | [**read\_init\_params**](classNyx.md#function-read-init-params-2-2) () <br>_Read initialization-related inputs._  |
|  void | [**read\_particle\_params**](classNyx.md#function-read-particle-params-1-2) () <br>_Read particle-related inputs._  |
|  void | [**read\_particle\_params**](classNyx.md#function-read-particle-params-2-2) () <br>_Read particle-related inputs._  |
|  void | [**set\_simd\_width**](classNyx.md#function-set-simd-width-1-2) (const int simd\_width) <br> |
|  void | [**set\_simd\_width**](classNyx.md#function-set-simd-width-1-2) (const int simd\_width) <br> |
|  amrex::Vector&lt; [**NyxParticleContainerBase**](classNyxParticleContainerBase.md) \* &gt; & | [**theActiveParticles**](classNyx.md#function-theactiveparticles-1-2) () <br> |
|  amrex::Vector&lt; [**NyxParticleContainerBase**](classNyxParticleContainerBase.md) \* &gt; & | [**theActiveParticles**](classNyx.md#function-theactiveparticles-2-2) () <br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**theDMPC**](classNyx.md#function-thedmpc-1-2) () <br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**theDMPC**](classNyx.md#function-thedmpc-2-2) () <br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**theGhostPC**](classNyx.md#function-theghostpc-1-2) () <br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**theGhostPC**](classNyx.md#function-theghostpc-2-2) () <br> |
|  amrex::Vector&lt; [**NyxParticleContainerBase**](classNyxParticleContainerBase.md) \* &gt; & | [**theGhostParticles**](classNyx.md#function-theghostparticles-1-2) () <br> |
|  amrex::Vector&lt; [**NyxParticleContainerBase**](classNyxParticleContainerBase.md) \* &gt; & | [**theGhostParticles**](classNyx.md#function-theghostparticles-2-2) () <br> |
|  [**StellarParticleContainer**](classNyxParticleContainer.md) \* | [**theGhostSPC**](classNyx.md#function-theghostspc-1-2) () <br> |
|  [**StellarParticleContainer**](classNyxParticleContainer.md) \* | [**theGhostSPC**](classNyx.md#function-theghostspc-2-2) () <br> |
|  [**StellarParticleContainer**](classNyxParticleContainer.md) \* | [**theSPC**](classNyx.md#function-thespc-1-2) () <br> |
|  [**StellarParticleContainer**](classNyxParticleContainer.md) \* | [**theSPC**](classNyx.md#function-thespc-2-2) () <br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**theVirtPC**](classNyx.md#function-thevirtpc-1-2) () <br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**theVirtPC**](classNyx.md#function-thevirtpc-2-2) () <br> |
|  [**StellarParticleContainer**](classNyxParticleContainer.md) \* | [**theVirtSPC**](classNyx.md#function-thevirtspc-1-2) () <br> |
|  [**StellarParticleContainer**](classNyxParticleContainer.md) \* | [**theVirtSPC**](classNyx.md#function-thevirtspc-2-2) () <br> |
|  amrex::Vector&lt; [**NyxParticleContainerBase**](classNyxParticleContainerBase.md) \* &gt; & | [**theVirtualParticles**](classNyx.md#function-thevirtualparticles-1-2) () <br> |
|  amrex::Vector&lt; [**NyxParticleContainerBase**](classNyxParticleContainerBase.md) \* &gt; & | [**theVirtualParticles**](classNyx.md#function-thevirtualparticles-2-2) () <br> |
|  void | [**variable\_cleanup**](classNyx.md#function-variable-cleanup-1-2) () <br>_Cleanup data descriptors at end of run._  |
|  void | [**variable\_cleanup**](classNyx.md#function-variable-cleanup-2-2) () <br>_Cleanup data descriptors at end of run._  |
|  void | [**variable\_setup**](classNyx.md#function-variable-setup-1-2) () <br>_Define data descriptors._  |
|  void | [**variable\_setup**](classNyx.md#function-variable-setup-2-2) () <br>_Define data descriptors._  |

## Public Static Functions inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

| Type | Name |
| ---: | :--- |
|  void | [**FillPatch**](classamrex_1_1AmrLevel.md#function-fillpatch) ([**AmrLevel**](classamrex_1_1AmrLevel.md) & amrlevel, MultiFab & leveldata, int boxGrow, Real time, int index, int scomp, int ncomp, int dcomp=0) <br> |
|  void | [**FillPatchAdd**](classamrex_1_1AmrLevel.md#function-fillpatchadd) ([**AmrLevel**](classamrex_1_1AmrLevel.md) & amrlevel, MultiFab & leveldata, int boxGrow, Real time, int index, int scomp, int ncomp, int dcomp=0) <br> |
|  void | [**FlushFPICache**](classamrex_1_1AmrLevel.md#function-flushfpicache) () <br> |
|  DeriveList & | [**get\_derive\_lst**](classamrex_1_1AmrLevel.md#function-get-derive-lst) () noexcept<br>_Returns list of derived variables._  |
|  const DescriptorList & | [**get\_desc\_lst**](classamrex_1_1AmrLevel.md#function-get-desc-lst) () noexcept<br>_Returns list of Descriptors._  |
|  bool | [**isStateVariable**](classamrex_1_1AmrLevel.md#function-isstatevariable) (const std::string & name, int & state\_indx, int & ncomp) <br>_Is name a state variable?_  |

## Public Static Functions inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

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
|  bool | [**FillPatchedOldState\_ok**](classNyx.md#variable-fillpatchedoldstate-ok)  <br>_These are the value of "z" at which to perform analysis._  |
|  amrex::FluxRegister \* | [**flux\_reg**](classNyx.md#variable-flux-reg)  <br> |

## Protected Attributes inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

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

## Protected Attributes inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

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
|  int | [**NUM\_GROW**](classNyx.md#variable-num-grow)   = = -1<br> |
|  amrex::IntVect | [**Nrep**](classNyx.md#variable-nrep)  <br> |
|  int | [**add\_ext\_src**](classNyx.md#variable-add-ext-src)   = = 0<br>_if true, define an additional source term_  |
|  int | [**allow\_untagging**](classNyx.md#variable-allow-untagging)   = = 0<br> |
|  amrex::Vector&lt; amrex::Real &gt; | [**analysis\_z\_values**](classNyx.md#variable-analysis-z-values)  <br>_These are the value of "z" at which to dump plotfiles._  |
|  amrex::Real | [**average\_dm\_density**](classNyx.md#variable-average-dm-density)   = = 0<br> |
|  amrex::Real | [**average\_gas\_density**](classNyx.md#variable-average-gas-density)   = = 0<br> |
|  amrex::Real | [**average\_neutr\_density**](classNyx.md#variable-average-neutr-density)   = = 0<br> |
|  amrex::Real | [**average\_total\_density**](classNyx.md#variable-average-total-density)   = = 0<br> |
|  amrex::Real | [**cfl**](classNyx.md#variable-cfl)   = = 0.8<br> |
|  amrex::Real | [**change\_max**](classNyx.md#variable-change-max)   = = 1.1<br> |
|  int | [**corner\_coupling**](classNyx.md#variable-corner-coupling)   = = 1<br> |
|  bool | [**do\_dm\_particles**](classNyx.md#variable-do-dm-particles)   = = false<br> |
|  int | [**do\_forcing**](classNyx.md#variable-do-forcing)   = =  0<br>_permits forcing to be switched on and off_  |
|  int | [**do\_grav**](classNyx.md#variable-do-grav)   = =  0<br>_permits gravity calculation to be turned on and off_  |
|  int | [**do\_hydro**](classNyx.md#variable-do-hydro)   = = -1<br>_permits hydro to be turned on and off for running pure rad problems:_  |
|  int | [**do\_reflux**](classNyx.md#variable-do-reflux)   = = 1<br> |
|  int | [**do\_special\_tagging**](classNyx.md#variable-do-special-tagging)   = = 0<br> |
|  bool | [**dump\_old**](classNyx.md#variable-dump-old)   = = false<br> |
|  amrex::ErrorList | [**err\_list**](classNyx.md#variable-err-list)  <br> |
|  amrex::Real | [**gamma**](classNyx.md#variable-gamma)   = =  0<br> |
|  amrex::Real | [**h\_species**](classNyx.md#variable-h-species)   = = 0.0<br> |
|  amrex::Real | [**he\_species**](classNyx.md#variable-he-species)   = = 0.0<br> |
|  int | [**heat\_cool\_type**](classNyx.md#variable-heat-cool-type)   = = 0<br>_specifies the heating/cooling source term_  |
|  int | [**inhomo\_grid**](classNyx.md#variable-inhomo-grid)   = = -1<br> |
|  int | [**inhomo\_reion**](classNyx.md#variable-inhomo-reion)   = = 0<br>_specifies inhomogeneous reionization type_  |
|  std::string | [**inhomo\_zhi\_file**](classNyx.md#variable-inhomo-zhi-file)   = = ""<br> |
|  amrex::Real | [**init\_shrink**](classNyx.md#variable-init-shrink)   = = 1.0<br> |
|  int | [**normalize\_species**](classNyx.md#variable-normalize-species)   = = 0<br> |
|  int | [**nsteps\_from\_plotfile**](classNyx.md#variable-nsteps-from-plotfile)   = = -1<br> |
|  int | [**num\_particle\_ghosts**](classNyx.md#variable-num-particle-ghosts)   = = 1<br> |
|  std::string | [**particle\_init\_type**](classNyx.md#variable-particle-init-type)   = = ""<br> |
|  long | [**particle\_initrandom\_count**](classNyx.md#variable-particle-initrandom-count)  <br> |
|  long | [**particle\_initrandom\_count\_per\_box**](classNyx.md#variable-particle-initrandom-count-per-box)  <br> |
|  int | [**particle\_initrandom\_iseed**](classNyx.md#variable-particle-initrandom-iseed)  <br> |
|  amrex::Real | [**particle\_initrandom\_mass**](classNyx.md#variable-particle-initrandom-mass)  <br> |
|  bool | [**particle\_initrandom\_serialize**](classNyx.md#variable-particle-initrandom-serialize)   = = false<br> |
|  std::string | [**particle\_move\_type**](classNyx.md#variable-particle-move-type)   = = ""<br> |
|  int | [**particle\_skip\_factor**](classNyx.md#variable-particle-skip-factor)   = = 1<br> |
|  amrex::BCRec | [**phys\_bc**](classNyx.md#variable-phys-bc)  <br> |
|  amrex::Vector&lt; amrex::Real &gt; | [**plot\_z\_values**](classNyx.md#variable-plot-z-values)  <br>_how many times the initial conditions are replicated in each direction_  |
|  int | [**ppm\_flatten\_before\_integrals**](classNyx.md#variable-ppm-flatten-before-integrals)   = = 0<br> |
|  int | [**ppm\_reference**](classNyx.md#variable-ppm-reference)   = = 1<br> |
|  int | [**ppm\_type**](classNyx.md#variable-ppm-type)   = = 1<br> |
|  amrex::Real | [**previousCPUTimeUsed**](classNyx.md#variable-previouscputimeused)   = = 0.0<br> |
|  amrex::Real | [**small\_dens**](classNyx.md#variable-small-dens)   = = -1.e200<br> |
|  amrex::Real | [**small\_temp**](classNyx.md#variable-small-temp)   = = -1.e200<br> |
|  amrex::Real | [**startCPUTime**](classNyx.md#variable-startcputime)   = = 0.0<br> |
|  int | [**strang\_split**](classNyx.md#variable-strang-split)   = = 1<br>_if true , incorporate the source term through Strang-splitting_  |
|  int | [**use\_colglaz**](classNyx.md#variable-use-colglaz)   = = 0<br> |
|  int | [**use\_const\_species**](classNyx.md#variable-use-const-species)   = = 0<br> |
|  int | [**use\_exact\_gravity**](classNyx.md#variable-use-exact-gravity)   = = 0<br> |
|  int | [**use\_flattening**](classNyx.md#variable-use-flattening)   = = 1<br> |
|  int | [**verbose**](classNyx.md#variable-verbose)   = = 0<br> |
|  int | [**version\_2**](classNyx.md#variable-version-2)   = = 0<br> |

## Protected Static Attributes inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

| Type | Name |
| ---: | :--- |
|  DeriveList | [**derive\_lst**](classamrex_1_1AmrLevel.md#variable-derive-lst)  <br> |
|  DescriptorList | [**desc\_lst**](classamrex_1_1AmrLevel.md#variable-desc-lst)  <br> |

## Protected Static Attributes inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

| Type | Name |
| ---: | :--- |
|  DeriveList | [**derive\_lst**](classamrex_1_1AmrLevel.md#variable-derive-lst)  <br> |
|  DescriptorList | [**desc\_lst**](classamrex_1_1AmrLevel.md#variable-desc-lst)  <br> |

## Protected Functions

| Type | Name |
| ---: | :--- |
|  void | [**average\_down**](classNyx.md#function-average-down-1-4) () <br> |
|  void | [**average\_down**](classNyx.md#function-average-down-2-4) (int state\_indx) <br> |
|  void | [**average\_down**](classNyx.md#function-average-down-1-4) () <br> |
|  void | [**average\_down**](classNyx.md#function-average-down-2-4) (int state\_indx) <br> |
|  void | [**build\_metrics**](classNyx.md#function-build-metrics-1-2) () <br> |
|  void | [**build\_metrics**](classNyx.md#function-build-metrics-1-2) () <br> |
|  void | [**compute\_average\_density**](classNyx.md#function-compute-average-density-1-2) () <br> |
|  void | [**compute\_average\_density**](classNyx.md#function-compute-average-density-1-2) () <br> |
|  void | [**compute\_average\_species**](classNyx.md#function-compute-average-species-1-2) (int nspec, int naux, amrex::Vector&lt; amrex::Real &gt; & average\_species) <br> |
|  void | [**compute\_average\_species**](classNyx.md#function-compute-average-species-1-2) (int nspec, int naux, amrex::Vector&lt; amrex::Real &gt; & average\_species) <br> |
|  void | [**compute\_average\_temperature**](classNyx.md#function-compute-average-temperature-1-2) (amrex::Real & average\_temperature) <br> |
|  void | [**compute\_average\_temperature**](classNyx.md#function-compute-average-temperature-1-2) (amrex::Real & average\_temperature) <br> |
|  void | [**enforce\_consistent\_e**](classNyx.md#function-enforce-consistent-e-1-2) (amrex::MultiFab & S) <br> |
|  void | [**enforce\_consistent\_e**](classNyx.md#function-enforce-consistent-e-1-2) (amrex::MultiFab & S) <br> |
|  void | [**enforce\_nonnegative\_species**](classNyx.md#function-enforce-nonnegative-species-1-2) (amrex::MultiFab & S\_new) <br> |
|  void | [**enforce\_nonnegative\_species**](classNyx.md#function-enforce-nonnegative-species-1-2) (amrex::MultiFab & S\_new) <br> |
|  amrex::FluxRegister & | [**get\_flux\_reg**](classNyx.md#function-get-flux-reg-1-4) () <br> |
|  amrex::FluxRegister & | [**get\_flux\_reg**](classNyx.md#function-get-flux-reg-2-4) (int lev) <br> |
|  amrex::FluxRegister & | [**get\_flux\_reg**](classNyx.md#function-get-flux-reg-3-4) () <br> |
|  amrex::FluxRegister & | [**get\_flux\_reg**](classNyx.md#function-get-flux-reg-4-4) (int lev) <br> |
|  [**Nyx**](classNyx.md) & | [**get\_level**](classNyx.md#function-get-level-1-2) (int lev) <br> |
|  [**Nyx**](classNyx.md) & | [**get\_level**](classNyx.md#function-get-level-2-2) (int lev) <br> |
|  void | [**reflux**](classNyx.md#function-reflux-1-2) () <br> |
|  void | [**reflux**](classNyx.md#function-reflux-1-2) () <br> |
|  std::string | [**retrieveDM**](classNyx.md#function-retrievedm-1-2) () <br> |
|  std::string | [**retrieveDM**](classNyx.md#function-retrievedm-1-2) () <br> |
|  void | [**set\_small\_values**](classNyx.md#function-set-small-values-1-2) () <br> |
|  void | [**set\_small\_values**](classNyx.md#function-set-small-values-1-2) () <br> |
| virtual void | [**sum\_integrated\_quantities**](classNyx.md#function-sum-integrated-quantities-1-2) () <br> |
| virtual void | [**sum\_integrated\_quantities**](classNyx.md#function-sum-integrated-quantities-2-2) () <br> |
|  void | [**write\_info**](classNyx.md#function-write-info-1-2) () <br> |
|  void | [**write\_info**](classNyx.md#function-write-info-1-2) () <br> |

## Protected Functions inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

| Type | Name |
| ---: | :--- |
|   | [**AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel-1-3) () noexcept<br>_The constructors_  _for derived classes._ |
|   | [**AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel-2-3) ([**Amr**](classamrex_1_1Amr.md) & papa, int lev, const Geometry & level\_geom, const BoxArray & bl, const DistributionMapping & dm, Real time) <br> |
|   | [**AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel-3-3) (const [**AmrLevel**](classamrex_1_1AmrLevel.md) &) = delete<br> |
|  void | [**finishConstructor**](classamrex_1_1AmrLevel.md#function-finishconstructor) () <br>_Common code used by all constructors._  |
|  [**AmrLevel**](classamrex_1_1AmrLevel.md) & | [**operator=**](classamrex_1_1AmrLevel.md#function-operator) (const [**AmrLevel**](classamrex_1_1AmrLevel.md) &) = delete<br> |

## Protected Functions inherited from [amrex::AmrLevel](classamrex_1_1AmrLevel.md)

| Type | Name |
| ---: | :--- |
|   | [**AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel-1-3) () noexcept<br>_The constructors_  _for derived classes._ |
|   | [**AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel-2-3) ([**Amr**](classamrex_1_1Amr.md) & papa, int lev, const Geometry & level\_geom, const BoxArray & bl, const DistributionMapping & dm, Real time) <br> |
|   | [**AmrLevel**](classamrex_1_1AmrLevel.md#function-amrlevel-3-3) (const [**AmrLevel**](classamrex_1_1AmrLevel.md) &) = delete<br> |
|  void | [**finishConstructor**](classamrex_1_1AmrLevel.md#function-finishconstructor) () <br>_Common code used by all constructors._  |
|  [**AmrLevel**](classamrex_1_1AmrLevel.md) & | [**operator=**](classamrex_1_1AmrLevel.md#function-operator) (const [**AmrLevel**](classamrex_1_1AmrLevel.md) &) = delete<br> |

## Protected Static Functions

| Type | Name |
| ---: | :--- |
|  amrex::Real | [**getCPUTime**](classNyx.md#function-getcputime-1-2) () <br> |
|  amrex::Real | [**getCPUTime**](classNyx.md#function-getcputime-2-2) () <br> |
|  void | [**network\_init**](classNyx.md#function-network-init-1-2) () <br> |
|  void | [**network\_init**](classNyx.md#function-network-init-2-2) () <br> |
|  void | [**read\_params**](classNyx.md#function-read-params-1-2) () <br> |
|  void | [**read\_params**](classNyx.md#function-read-params-2-2) () <br> |



## Public Attributes Documentation


### <a href="#variable-fine-mask" id="variable-fine-mask">variable fine\_mask </a>


```cpp
amrex::MultiFab * Nyx::fine_mask;
```


## Public Static Attributes Documentation


### <a href="#variable-density" id="variable-density">variable Density </a>


```cpp
static int Nyx::Density;
```



### <a href="#variable-eden" id="variable-eden">variable Eden </a>


```cpp
static int Nyx::Eden;
```



### <a href="#variable-eint" id="variable-eint">variable Eint </a>


```cpp
static int Nyx::Eint;
```



### <a href="#variable-firstadv" id="variable-firstadv">variable FirstAdv </a>


```cpp
static int Nyx::FirstAdv;
```



### <a href="#variable-firstaux" id="variable-firstaux">variable FirstAux </a>


```cpp
static int Nyx::FirstAux;
```



### <a href="#variable-firstspec" id="variable-firstspec">variable FirstSpec </a>


```cpp
static int Nyx::FirstSpec;
```



### <a href="#variable-num-state" id="variable-num-state">variable NUM\_STATE </a>


```cpp
static int Nyx::NUM_STATE;
```



### <a href="#variable-ne-comp" id="variable-ne-comp">variable Ne\_comp </a>


```cpp
static int Nyx::Ne_comp;
```



### <a href="#variable-numadv" id="variable-numadv">variable NumAdv </a>


```cpp
static int Nyx::NumAdv;
```



### <a href="#variable-numaux" id="variable-numaux">variable NumAux </a>


```cpp
static int Nyx::NumAux;
```



### <a href="#variable-numspec" id="variable-numspec">variable NumSpec </a>


```cpp
static int Nyx::NumSpec;
```



### <a href="#variable-temp-comp" id="variable-temp-comp">variable Temp\_comp </a>


```cpp
static int Nyx::Temp_comp;
```



### <a href="#variable-xmom" id="variable-xmom">variable Xmom </a>


```cpp
static int Nyx::Xmom;
```



### <a href="#variable-ymom" id="variable-ymom">variable Ymom </a>


```cpp
static int Nyx::Ymom;
```



### <a href="#variable-zhi-comp" id="variable-zhi-comp">variable Zhi\_comp </a>


```cpp
static int Nyx::Zhi_comp;
```



### <a href="#variable-zmom" id="variable-zmom">variable Zmom </a>


```cpp
static int Nyx::Zmom;
```



### <a href="#variable-absolute-max-change-a" id="variable-absolute-max-change-a">variable absolute\_max\_change\_a </a>


```cpp
static amrex::Real Nyx::absolute_max_change_a;
```



### <a href="#variable-comoving-omb" id="variable-comoving-omb">variable comoving\_OmB </a>


```cpp
static amrex::Real Nyx::comoving_OmB;
```



### <a href="#variable-comoving-omm" id="variable-comoving-omm">variable comoving\_OmM </a>


```cpp
static amrex::Real Nyx::comoving_OmM;
```



### <a href="#variable-comoving-h" id="variable-comoving-h">variable comoving\_h </a>


```cpp
static amrex::Real Nyx::comoving_h;
```



### <a href="#variable-dt-binpow" id="variable-dt-binpow">variable dt\_binpow </a>


```cpp
static amrex::Real Nyx::dt_binpow;
```



### <a href="#variable-final-a" id="variable-final-a">variable final\_a </a>


```cpp
static amrex::Real Nyx::final_a;
```



### <a href="#variable-final-time" id="variable-final-time">variable final\_time </a>


```cpp
static amrex::Real Nyx::final_time;
```



### <a href="#variable-final-z" id="variable-final-z">variable final\_z </a>


```cpp
static amrex::Real Nyx::final_z;
```



### <a href="#variable-init-with-sph-particles" id="variable-init-with-sph-particles">variable init\_with\_sph\_particles </a>


```cpp
static int Nyx::init_with_sph_particles;
```



### <a href="#variable-initial-time" id="variable-initial-time">variable initial\_time </a>


```cpp
static amrex::Real Nyx::initial_time;
```



### <a href="#variable-initial-z" id="variable-initial-z">variable initial\_z </a>


```cpp
static amrex::Real Nyx::initial_z;
```



### <a href="#variable-new-a" id="variable-new-a">variable new\_a </a>


```cpp
static amrex::Real Nyx::new_a;
```



### <a href="#variable-new-a-time" id="variable-new-a-time">variable new\_a\_time </a>


```cpp
static amrex::Real Nyx::new_a_time;
```



### <a href="#variable-old-a" id="variable-old-a">variable old\_a </a>


```cpp
static amrex::Real Nyx::old_a;
```



### <a href="#variable-old-a-time" id="variable-old-a-time">variable old\_a\_time </a>


```cpp
static amrex::Real Nyx::old_a_time;
```



### <a href="#variable-particle-cfl" id="variable-particle-cfl">variable particle\_cfl </a>


```cpp
Real Nyx::particle_cfl;
```



### <a href="#variable-particle-plotfile-format" id="variable-particle-plotfile-format">variable particle\_plotfile\_format </a>


```cpp
static std::string Nyx::particle_plotfile_format;
```



### <a href="#variable-particle-verbose" id="variable-particle-verbose">variable particle\_verbose </a>


```cpp
int Nyx::particle_verbose;
```



### <a href="#variable-print-fortran-warnings" id="variable-print-fortran-warnings">variable print\_fortran\_warnings </a>


```cpp
static int Nyx::print_fortran_warnings;
```



### <a href="#variable-relative-max-change-a" id="variable-relative-max-change-a">variable relative\_max\_change\_a </a>


```cpp
static amrex::Real Nyx::relative_max_change_a;
```



### <a href="#variable-strict-subcycling" id="variable-strict-subcycling">variable strict\_subcycling </a>


```cpp
static int Nyx::strict_subcycling;
```



### <a href="#variable-write-coarsened-particles" id="variable-write-coarsened-particles">variable write\_coarsened\_particles </a>


```cpp
int Nyx::write_coarsened_particles;
```


Shall we write an ascii file with only the particles in every other cell  this is a cheap way of creating a "coarser" particle file 


        

### <a href="#variable-write-parameters-in-plotfile" id="variable-write-parameters-in-plotfile">variable write\_parameters\_in\_plotfile </a>


```cpp
static int Nyx::write_parameters_in_plotfile;
```



### <a href="#variable-write-particle-density-at-init" id="variable-write-particle-density-at-init">variable write\_particle\_density\_at\_init </a>


```cpp
int Nyx::write_particle_density_at_init;
```


Shall we write the initial single-level particle density into a multifab called "ParticleDensity"? 


        
## Public Functions Documentation


### <a href="#function-createleveldirectory-1-2" id="function-createleveldirectory-1-2">function CreateLevelDirectory [1/2]</a>


```cpp
virtual void Nyx::CreateLevelDirectory (
    const std::string & dir
) 
```


Implements [*amrex::AmrLevel::CreateLevelDirectory*](classamrex_1_1AmrLevel.md#function-createleveldirectory)


### <a href="#function-createleveldirectory-2-2" id="function-createleveldirectory-2-2">function CreateLevelDirectory [2/2]</a>


```cpp
virtual void Nyx::CreateLevelDirectory (
    const std::string & dir
) 
```


Implements [*amrex::AmrLevel::CreateLevelDirectory*](classamrex_1_1AmrLevel.md#function-createleveldirectory)


### <a href="#function-leveldirectorynames-1-2" id="function-leveldirectorynames-1-2">function LevelDirectoryNames [1/2]</a>


```cpp
void Nyx::LevelDirectoryNames (
    const std::string & dir,
    const std::string & secondDir,
    std::string & LevelDir,
    std::string & FullPath
) 
```



### <a href="#function-leveldirectorynames-1-2" id="function-leveldirectorynames-1-2">function LevelDirectoryNames [1/2]</a>


```cpp
void Nyx::LevelDirectoryNames (
    const std::string & dir,
    const std::string & secondDir,
    std::string & LevelDir,
    std::string & FullPath
) 
```



### <a href="#function-lya-statistics-1-2" id="function-lya-statistics-1-2">function Lya\_statistics [1/2]</a>


```cpp
void Nyx::Lya_statistics () 
```



### <a href="#function-lya-statistics-1-2" id="function-lya-statistics-1-2">function Lya\_statistics [1/2]</a>


```cpp
void Nyx::Lya_statistics () 
```



### <a href="#function-nyx-1-4" id="function-nyx-1-4">function Nyx [1/4]</a>


```cpp
Nyx::Nyx () 
```



### <a href="#function-nyx-2-4" id="function-nyx-2-4">function Nyx [2/4]</a>


```cpp
Nyx::Nyx (
    amrex::Amr & papa,
    int lev,
    const amrex::Geometry & level_geom,
    const amrex::BoxArray & bl,
    const amrex::DistributionMapping & dm,
    amrex::Real time
) 
```



### <a href="#function-nyx-1-4" id="function-nyx-1-4">function Nyx [1/4]</a>


```cpp
Nyx::Nyx () 
```



### <a href="#function-nyx-2-4" id="function-nyx-2-4">function Nyx [2/4]</a>


```cpp
Nyx::Nyx (
    amrex::Amr & papa,
    int lev,
    const amrex::Geometry & level_geom,
    const amrex::BoxArray & bl,
    const amrex::DistributionMapping & dm,
    amrex::Real time
) 
```



### <a href="#function-readplotfile-1-2" id="function-readplotfile-1-2">function ReadPlotFile [1/2]</a>


```cpp
void Nyx::ReadPlotFile (
    bool first,
    const std::string & plot_file_name,
    bool & rhoe_infile
) 
```



### <a href="#function-readplotfile-1-2" id="function-readplotfile-1-2">function ReadPlotFile [1/2]</a>


```cpp
void Nyx::ReadPlotFile (
    bool first,
    const std::string & plot_file_name,
    bool & rhoe_infile
) 
```



### <a href="#function-advance-1-2" id="function-advance-1-2">function advance [1/2]</a>


```cpp
virtual amrex::Real Nyx::advance (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-advance-1-2" id="function-advance-1-2">function advance [1/2]</a>


```cpp
virtual amrex::Real Nyx::advance (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-advance-hydro-1-2" id="function-advance-hydro-1-2">function advance\_hydro [1/2]</a>


```cpp
amrex::Real Nyx::advance_hydro (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-advance-hydro-1-2" id="function-advance-hydro-1-2">function advance\_hydro [1/2]</a>


```cpp
amrex::Real Nyx::advance_hydro (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-advance-hydro-plus-particles-1-2" id="function-advance-hydro-plus-particles-1-2">function advance\_hydro\_plus\_particles [1/2]</a>


```cpp
amrex::Real Nyx::advance_hydro_plus_particles (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-advance-hydro-plus-particles-1-2" id="function-advance-hydro-plus-particles-1-2">function advance\_hydro\_plus\_particles [1/2]</a>


```cpp
amrex::Real Nyx::advance_hydro_plus_particles (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-advance-no-hydro-1-2" id="function-advance-no-hydro-1-2">function advance\_no\_hydro [1/2]</a>


```cpp
amrex::Real Nyx::advance_no_hydro (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-advance-no-hydro-1-2" id="function-advance-no-hydro-1-2">function advance\_no\_hydro [1/2]</a>


```cpp
amrex::Real Nyx::advance_no_hydro (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-advance-particles-only-1-2" id="function-advance-particles-only-1-2">function advance\_particles\_only [1/2]</a>


```cpp
amrex::Real Nyx::advance_particles_only (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-advance-particles-only-1-2" id="function-advance-particles-only-1-2">function advance\_particles\_only [1/2]</a>


```cpp
amrex::Real Nyx::advance_particles_only (
    amrex::Real time,
    amrex::Real dt,
    int iteration,
    int ncycle
) 
```



### <a href="#function-analysis-z-est-time-step-1-2" id="function-analysis-z-est-time-step-1-2">function analysis\_z\_est\_time\_step [1/2]</a>


```cpp
void Nyx::analysis_z_est_time_step (
    amrex::Real & est_dt,
    bool & dt_changed
) 
```



### <a href="#function-analysis-z-est-time-step-1-2" id="function-analysis-z-est-time-step-1-2">function analysis\_z\_est\_time\_step [1/2]</a>


```cpp
void Nyx::analysis_z_est_time_step (
    amrex::Real & est_dt,
    bool & dt_changed
) 
```



### <a href="#function-build-fine-mask-1-2" id="function-build-fine-mask-1-2">function build\_fine\_mask [1/2]</a>


```cpp
amrex::MultiFab * Nyx::build_fine_mask () 
```



### <a href="#function-build-fine-mask-2-2" id="function-build-fine-mask-2-2">function build\_fine\_mask [2/2]</a>


```cpp
amrex::MultiFab * Nyx::build_fine_mask () 
```



### <a href="#function-checkpoint-1-2" id="function-checkpoint-1-2">function checkPoint [1/2]</a>


```cpp
virtual void Nyx::checkPoint (
    const std::string & dir,
    std::ostream & os,
    amrex::VisMF::How how,
    bool dump_old
) 
```



### <a href="#function-checkpoint-1-2" id="function-checkpoint-1-2">function checkPoint [1/2]</a>


```cpp
virtual void Nyx::checkPoint (
    const std::string & dir,
    std::ostream & os,
    amrex::VisMF::How how,
    bool dump_old
) 
```



### <a href="#function-checkpointpost-1-2" id="function-checkpointpost-1-2">function checkPointPost [1/2]</a>


```cpp
virtual void Nyx::checkPointPost (
    const std::string & dir,
    std::ostream & os
) 
```


Implements [*amrex::AmrLevel::checkPointPost*](classamrex_1_1AmrLevel.md#function-checkpointpost)


### <a href="#function-checkpointpost-2-2" id="function-checkpointpost-2-2">function checkPointPost [2/2]</a>


```cpp
virtual void Nyx::checkPointPost (
    const std::string & dir,
    std::ostream & os
) 
```


Implements [*amrex::AmrLevel::checkPointPost*](classamrex_1_1AmrLevel.md#function-checkpointpost)


### <a href="#function-checkpointpre-1-2" id="function-checkpointpre-1-2">function checkPointPre [1/2]</a>


```cpp
virtual void Nyx::checkPointPre (
    const std::string & dir,
    std::ostream & os
) 
```


Implements [*amrex::AmrLevel::checkPointPre*](classamrex_1_1AmrLevel.md#function-checkpointpre)


### <a href="#function-checkpointpre-2-2" id="function-checkpointpre-2-2">function checkPointPre [2/2]</a>


```cpp
virtual void Nyx::checkPointPre (
    const std::string & dir,
    std::ostream & os
) 
```


Implements [*amrex::AmrLevel::checkPointPre*](classamrex_1_1AmrLevel.md#function-checkpointpre)


### <a href="#function-comoving-a-post-restart-1-2" id="function-comoving-a-post-restart-1-2">function comoving\_a\_post\_restart [1/2]</a>


```cpp
void Nyx::comoving_a_post_restart (
    const std::string & restart_file
) 
```



### <a href="#function-comoving-a-post-restart-1-2" id="function-comoving-a-post-restart-1-2">function comoving\_a\_post\_restart [1/2]</a>


```cpp
void Nyx::comoving_a_post_restart (
    const std::string & restart_file
) 
```



### <a href="#function-comoving-est-time-step-1-2" id="function-comoving-est-time-step-1-2">function comoving\_est\_time\_step [1/2]</a>


```cpp
void Nyx::comoving_est_time_step (
    amrex::Real & cur_time,
    amrex::Real & est_dt
) 
```



### <a href="#function-comoving-est-time-step-1-2" id="function-comoving-est-time-step-1-2">function comoving\_est\_time\_step [1/2]</a>


```cpp
void Nyx::comoving_est_time_step (
    amrex::Real & cur_time,
    amrex::Real & est_dt
) 
```



### <a href="#function-computeinitialdt-1-2" id="function-computeinitialdt-1-2">function computeInitialDt [1/2]</a>


```cpp
virtual void Nyx::computeInitialDt (
    int finest_level,
    int sub_cycle,
    amrex::Vector< int > & n_cycle,
    const amrex::Vector< amrex::IntVect > & ref_ratio,
    amrex::Vector< amrex::Real > & dt_level,
    amrex::Real stop_time
) 
```



### <a href="#function-computeinitialdt-1-2" id="function-computeinitialdt-1-2">function computeInitialDt [1/2]</a>


```cpp
virtual void Nyx::computeInitialDt (
    int finest_level,
    int sub_cycle,
    amrex::Vector< int > & n_cycle,
    const amrex::Vector< amrex::IntVect > & ref_ratio,
    amrex::Vector< amrex::Real > & dt_level,
    amrex::Real stop_time
) 
```



### <a href="#function-computenewdt-1-2" id="function-computenewdt-1-2">function computeNewDt [1/2]</a>


```cpp
virtual void Nyx::computeNewDt (
    int finest_level,
    int sub_cycle,
    amrex::Vector< int > & n_cycle,
    const amrex::Vector< amrex::IntVect > & ref_ratio,
    amrex::Vector< amrex::Real > & dt_min,
    amrex::Vector< amrex::Real > & dt_level,
    amrex::Real stop_time,
    int post_regrid_flag
) 
```



### <a href="#function-computenewdt-1-2" id="function-computenewdt-1-2">function computeNewDt [1/2]</a>


```cpp
virtual void Nyx::computeNewDt (
    int finest_level,
    int sub_cycle,
    amrex::Vector< int > & n_cycle,
    const amrex::Vector< amrex::IntVect > & ref_ratio,
    amrex::Vector< amrex::Real > & dt_min,
    amrex::Vector< amrex::Real > & dt_level,
    amrex::Real stop_time,
    int post_regrid_flag
) 
```



### <a href="#function-compute-gas-fractions-1-2" id="function-compute-gas-fractions-1-2">function compute\_gas\_fractions [1/2]</a>


```cpp
void Nyx::compute_gas_fractions (
    amrex::Real T_cut,
    amrex::Real rho_cut,
    amrex::Real & whim_mass_frac,
    amrex::Real & whim_vol_frac,
    amrex::Real & hh_mass_frac,
    amrex::Real & hh_vol_frac,
    amrex::Real & igm_mass_frac,
    amrex::Real & igm_vol_frac
) 
```



### <a href="#function-compute-gas-fractions-1-2" id="function-compute-gas-fractions-1-2">function compute\_gas\_fractions [1/2]</a>


```cpp
void Nyx::compute_gas_fractions (
    amrex::Real T_cut,
    amrex::Real rho_cut,
    amrex::Real & whim_mass_frac,
    amrex::Real & whim_vol_frac,
    amrex::Real & hh_mass_frac,
    amrex::Real & hh_vol_frac,
    amrex::Real & igm_mass_frac,
    amrex::Real & igm_vol_frac
) 
```



### <a href="#function-compute-hydro-sources-1-2" id="function-compute-hydro-sources-1-2">function compute\_hydro\_sources [1/2]</a>


```cpp
void Nyx::compute_hydro_sources (
    amrex::Real time,
    amrex::Real dt,
    amrex::Real a_old,
    amrex::Real a_new,
    amrex::MultiFab & S_border,
    amrex::MultiFab & D_border,
    amrex::MultiFab & ext_src_old,
    amrex::MultiFab & hydro_src,
    amrex::MultiFab & grav,
    amrex::MultiFab & divu_cc,
    bool init_flux_register,
    bool add_to_flux_register
) 
```



### <a href="#function-compute-hydro-sources-1-2" id="function-compute-hydro-sources-1-2">function compute\_hydro\_sources [1/2]</a>


```cpp
void Nyx::compute_hydro_sources (
    amrex::Real time,
    amrex::Real dt,
    amrex::Real a_old,
    amrex::Real a_new,
    amrex::MultiFab & S_border,
    amrex::MultiFab & D_border,
    amrex::MultiFab & ext_src_old,
    amrex::MultiFab & hydro_src,
    amrex::MultiFab & grav,
    amrex::MultiFab & divu_cc,
    bool init_flux_register,
    bool add_to_flux_register
) 
```



### <a href="#function-compute-new-temp-1-2" id="function-compute-new-temp-1-2">function compute\_new\_temp [1/2]</a>


```cpp
void Nyx::compute_new_temp (
    amrex::MultiFab & S_new,
    amrex::MultiFab & D_new
) 
```



### <a href="#function-compute-new-temp-1-2" id="function-compute-new-temp-1-2">function compute\_new\_temp [1/2]</a>


```cpp
void Nyx::compute_new_temp (
    amrex::MultiFab & S_new,
    amrex::MultiFab & D_new
) 
```



### <a href="#function-compute-rho-temp-1-2" id="function-compute-rho-temp-1-2">function compute\_rho\_temp [1/2]</a>


```cpp
void Nyx::compute_rho_temp (
    amrex::Real & rho_T_avg,
    amrex::Real & T_avg,
    amrex::Real & Tinv_avg,
    amrex::Real & T_meanrho
) 
```



### <a href="#function-compute-rho-temp-1-2" id="function-compute-rho-temp-1-2">function compute\_rho\_temp [1/2]</a>


```cpp
void Nyx::compute_rho_temp (
    amrex::Real & rho_T_avg,
    amrex::Real & T_avg,
    amrex::Real & Tinv_avg,
    amrex::Real & T_meanrho
) 
```



### <a href="#function-conserved-to-primitive-1-2" id="function-conserved-to-primitive-1-2">function conserved\_to\_primitive [1/2]</a>


```cpp
void Nyx::conserved_to_primitive (
    amrex::MultiFab & state
) 
```



### <a href="#function-conserved-to-primitive-1-2" id="function-conserved-to-primitive-1-2">function conserved\_to\_primitive [1/2]</a>


```cpp
void Nyx::conserved_to_primitive (
    amrex::MultiFab & state
) 
```



### <a href="#function-derive-1-4" id="function-derive-1-4">function derive [1/4]</a>


```cpp
std::unique_ptr< amrex::MultiFab > Nyx::derive (
    const std::string & name,
    amrex::Real time,
    int ngrow
) 
```


Returns a amrex::MultiFab containing the derived data for this level. The user is responsible for deleting this pointer when done with it. If `ngrow` &gt; 0 the amrex::MultiFab is built on the appropriately grown amrex::BoxArray. 


        

### <a href="#function-derive-2-4" id="function-derive-2-4">function derive [2/4]</a>


```cpp
void Nyx::derive (
    const std::string & name,
    amrex::Real time,
    amrex::MultiFab & mf,
    int dcomp
) 
```


This version of `derive()` fills the dcomp'th component of mf with the derived quantity. 


        

### <a href="#function-derive-1-4" id="function-derive-1-4">function derive [1/4]</a>


```cpp
std::unique_ptr< amrex::MultiFab > Nyx::derive (
    const std::string & name,
    amrex::Real time,
    int ngrow
) 
```


Returns a amrex::MultiFab containing the derived data for this level. The user is responsible for deleting this pointer when done with it. If `ngrow` &gt; 0 the amrex::MultiFab is built on the appropriately grown amrex::BoxArray. 


        

### <a href="#function-derive-2-4" id="function-derive-2-4">function derive [2/4]</a>


```cpp
void Nyx::derive (
    const std::string & name,
    amrex::Real time,
    amrex::MultiFab & mf,
    int dcomp
) 
```


This version of `derive()` fills the dcomp'th component of mf with the derived quantity. 


        

### <a href="#function-doanalysisnow-1-2" id="function-doanalysisnow-1-2">function doAnalysisNow [1/2]</a>


```cpp
bool Nyx::doAnalysisNow () 
```



### <a href="#function-doanalysisnow-1-2" id="function-doanalysisnow-1-2">function doAnalysisNow [1/2]</a>


```cpp
bool Nyx::doAnalysisNow () 
```



### <a href="#function-do-energy-diagnostics-1-2" id="function-do-energy-diagnostics-1-2">function do\_energy\_diagnostics [1/2]</a>


```cpp
void Nyx::do_energy_diagnostics () 
```



### <a href="#function-do-energy-diagnostics-1-2" id="function-do-energy-diagnostics-1-2">function do\_energy\_diagnostics [1/2]</a>


```cpp
void Nyx::do_energy_diagnostics () 
```



### <a href="#function-errorest-1-2" id="function-errorest-1-2">function errorEst [1/2]</a>


```cpp
virtual void Nyx::errorEst (
    amrex::TagBoxArray & tb,
    int clearval,
    int tagval,
    amrex::Real time,
    int n_error_buf=0,
    int ngrow=0
) 
```



### <a href="#function-errorest-1-2" id="function-errorest-1-2">function errorEst [1/2]</a>


```cpp
virtual void Nyx::errorEst (
    amrex::TagBoxArray & tb,
    int clearval,
    int tagval,
    amrex::Real time,
    int n_error_buf=0,
    int ngrow=0
) 
```



### <a href="#function-est-time-step-1-2" id="function-est-time-step-1-2">function est\_time\_step [1/2]</a>


```cpp
amrex::Real Nyx::est_time_step (
    amrex::Real dt_old
) 
```



### <a href="#function-est-time-step-1-2" id="function-est-time-step-1-2">function est\_time\_step [1/2]</a>


```cpp
amrex::Real Nyx::est_time_step (
    amrex::Real dt_old
) 
```



### <a href="#function-get-comoving-a-1-2" id="function-get-comoving-a-1-2">function get\_comoving\_a [1/2]</a>


```cpp
amrex::Real Nyx::get_comoving_a (
    amrex::Real time
) 
```



### <a href="#function-get-comoving-a-1-2" id="function-get-comoving-a-1-2">function get\_comoving\_a [1/2]</a>


```cpp
amrex::Real Nyx::get_comoving_a (
    amrex::Real time
) 
```



### <a href="#function-get-new-source-1-2" id="function-get-new-source-1-2">function get\_new\_source [1/2]</a>


```cpp
void Nyx::get_new_source (
    amrex::Real old_time,
    amrex::Real new_time,
    amrex::Real dt,
    amrex::MultiFab & Rhs
) 
```



### <a href="#function-get-new-source-1-2" id="function-get-new-source-1-2">function get\_new\_source [1/2]</a>


```cpp
void Nyx::get_new_source (
    amrex::Real old_time,
    amrex::Real new_time,
    amrex::Real dt,
    amrex::MultiFab & Rhs
) 
```



### <a href="#function-get-old-source-1-2" id="function-get-old-source-1-2">function get\_old\_source [1/2]</a>


```cpp
void Nyx::get_old_source (
    amrex::Real old_time,
    amrex::Real dt,
    amrex::MultiFab & Rhs
) 
```



### <a href="#function-get-old-source-1-2" id="function-get-old-source-1-2">function get\_old\_source [1/2]</a>


```cpp
void Nyx::get_old_source (
    amrex::Real old_time,
    amrex::Real dt,
    amrex::MultiFab & Rhs
) 
```



### <a href="#function-halo-accrete-1-2" id="function-halo-accrete-1-2">function halo\_accrete [1/2]</a>


```cpp
void Nyx::halo_accrete (
    amrex::Real dt
) 
```



### <a href="#function-halo-accrete-1-2" id="function-halo-accrete-1-2">function halo\_accrete [1/2]</a>


```cpp
void Nyx::halo_accrete (
    amrex::Real dt
) 
```



### <a href="#function-halo-find-1-2" id="function-halo-find-1-2">function halo\_find [1/2]</a>


```cpp
void Nyx::halo_find (
    amrex::Real dt
) 
```



### <a href="#function-halo-find-1-2" id="function-halo-find-1-2">function halo\_find [1/2]</a>


```cpp
void Nyx::halo_find (
    amrex::Real dt
) 
```



### <a href="#function-halo-merge-1-2" id="function-halo-merge-1-2">function halo\_merge [1/2]</a>


```cpp
void Nyx::halo_merge () 
```



### <a href="#function-halo-merge-1-2" id="function-halo-merge-1-2">function halo\_merge [1/2]</a>


```cpp
void Nyx::halo_merge () 
```



### <a href="#function-init-1-4" id="function-init-1-4">function init [1/4]</a>


```cpp
virtual void Nyx::init (
    amrex::AmrLevel & old
) 
```


Implements [*amrex::AmrLevel::init*](classamrex_1_1AmrLevel.md#function-init-1-2)


### <a href="#function-init-2-4" id="function-init-2-4">function init [2/4]</a>


```cpp
virtual void Nyx::init () 
```


Implements [*amrex::AmrLevel::init*](classamrex_1_1AmrLevel.md#function-init-2-2)


### <a href="#function-init-3-4" id="function-init-3-4">function init [3/4]</a>


```cpp
virtual void Nyx::init (
    amrex::AmrLevel & old
) 
```


Implements [*amrex::AmrLevel::init*](classamrex_1_1AmrLevel.md#function-init-1-2)


### <a href="#function-init-4-4" id="function-init-4-4">function init [4/4]</a>


```cpp
virtual void Nyx::init () 
```


Implements [*amrex::AmrLevel::init*](classamrex_1_1AmrLevel.md#function-init-2-2)


### <a href="#function-initdata-1-2" id="function-initdata-1-2">function initData [1/2]</a>


```cpp
virtual void Nyx::initData () 
```


Implements [*amrex::AmrLevel::initData*](classamrex_1_1AmrLevel.md#function-initdata)


### <a href="#function-initdata-2-2" id="function-initdata-2-2">function initData [2/2]</a>


```cpp
virtual void Nyx::initData () 
```


Implements [*amrex::AmrLevel::initData*](classamrex_1_1AmrLevel.md#function-initdata)


### <a href="#function-init-from-plotfile-1-2" id="function-init-from-plotfile-1-2">function init\_from\_plotfile [1/2]</a>


```cpp
void Nyx::init_from_plotfile () 
```



### <a href="#function-init-from-plotfile-1-2" id="function-init-from-plotfile-1-2">function init\_from\_plotfile [1/2]</a>


```cpp
void Nyx::init_from_plotfile () 
```



### <a href="#function-init-particles-1-2" id="function-init-particles-1-2">function init\_particles [1/2]</a>


```cpp
virtual void Nyx::init_particles () 
```



### <a href="#function-init-particles-2-2" id="function-init-particles-2-2">function init\_particles [2/2]</a>


```cpp
virtual void Nyx::init_particles () 
```



### <a href="#function-init-zhi-1-2" id="function-init-zhi-1-2">function init\_zhi [1/2]</a>


```cpp
void Nyx::init_zhi () 
```



### <a href="#function-init-zhi-1-2" id="function-init-zhi-1-2">function init\_zhi [1/2]</a>


```cpp
void Nyx::init_zhi () 
```



### <a href="#function-initcosmo-1-2" id="function-initcosmo-1-2">function initcosmo [1/2]</a>


```cpp
void Nyx::initcosmo () 
```



### <a href="#function-initcosmo-1-2" id="function-initcosmo-1-2">function initcosmo [1/2]</a>


```cpp
void Nyx::initcosmo () 
```



### <a href="#function-initial-time-step-1-2" id="function-initial-time-step-1-2">function initial\_time\_step [1/2]</a>


```cpp
amrex::Real Nyx::initial_time_step () 
```



### <a href="#function-initial-time-step-2-2" id="function-initial-time-step-2-2">function initial\_time\_step [2/2]</a>


```cpp
amrex::Real Nyx::initial_time_step () 
```



### <a href="#function-integrate-comoving-a-1-2" id="function-integrate-comoving-a-1-2">function integrate\_comoving\_a [1/2]</a>


```cpp
void Nyx::integrate_comoving_a (
    amrex::Real time,
    amrex::Real dt
) 
```



### <a href="#function-integrate-comoving-a-1-2" id="function-integrate-comoving-a-1-2">function integrate\_comoving\_a [1/2]</a>


```cpp
void Nyx::integrate_comoving_a (
    amrex::Real time,
    amrex::Real dt
) 
```



### <a href="#function-integrate-state-box-1-2" id="function-integrate-state-box-1-2">function integrate\_state\_box [1/2]</a>


```cpp
int Nyx::integrate_state_box (
    amrex::MultiFab & state,
    amrex::MultiFab & diag_eos,
    const amrex::Real & a,
    const amrex::Real & delta_time
) 
```



### <a href="#function-integrate-state-box-1-2" id="function-integrate-state-box-1-2">function integrate\_state\_box [1/2]</a>


```cpp
int Nyx::integrate_state_box (
    amrex::MultiFab & state,
    amrex::MultiFab & diag_eos,
    const amrex::Real & a,
    const amrex::Real & delta_time
) 
```



### <a href="#function-integrate-state-grownbox-1-2" id="function-integrate-state-grownbox-1-2">function integrate\_state\_grownbox [1/2]</a>


```cpp
int Nyx::integrate_state_grownbox (
    amrex::MultiFab & state,
    amrex::MultiFab & diag_eos,
    const amrex::Real & a,
    const amrex::Real & delta_time
) 
```



### <a href="#function-integrate-state-grownbox-1-2" id="function-integrate-state-grownbox-1-2">function integrate\_state\_grownbox [1/2]</a>


```cpp
int Nyx::integrate_state_grownbox (
    amrex::MultiFab & state,
    amrex::MultiFab & diag_eos,
    const amrex::Real & a,
    const amrex::Real & delta_time
) 
```



### <a href="#function-integrate-state-grownvec-1-2" id="function-integrate-state-grownvec-1-2">function integrate\_state\_grownvec [1/2]</a>


```cpp
int Nyx::integrate_state_grownvec (
    amrex::MultiFab & state,
    amrex::MultiFab & diag_eos,
    const amrex::Real & a,
    const amrex::Real & delta_time
) 
```



### <a href="#function-integrate-state-grownvec-1-2" id="function-integrate-state-grownvec-1-2">function integrate\_state\_grownvec [1/2]</a>


```cpp
int Nyx::integrate_state_grownvec (
    amrex::MultiFab & state,
    amrex::MultiFab & diag_eos,
    const amrex::Real & a,
    const amrex::Real & delta_time
) 
```



### <a href="#function-integrate-state-vec-1-2" id="function-integrate-state-vec-1-2">function integrate\_state\_vec [1/2]</a>


```cpp
int Nyx::integrate_state_vec (
    amrex::MultiFab & state,
    amrex::MultiFab & diag_eos,
    const amrex::Real & a,
    const amrex::Real & delta_time
) 
```



### <a href="#function-integrate-state-vec-1-2" id="function-integrate-state-vec-1-2">function integrate\_state\_vec [1/2]</a>


```cpp
int Nyx::integrate_state_vec (
    amrex::MultiFab & state,
    amrex::MultiFab & diag_eos,
    const amrex::Real & a,
    const amrex::Real & delta_time
) 
```



### <a href="#function-manual-tags-placement-1-2" id="function-manual-tags-placement-1-2">function manual\_tags\_placement [1/2]</a>


```cpp
virtual void Nyx::manual_tags_placement (
    amrex::TagBoxArray & tags,
    const amrex::Vector< amrex::IntVect > & bf_lev
) override
```


Called in grid\_places after other tagging routines to modify the list of tagged points 


        

### <a href="#function-manual-tags-placement-1-2" id="function-manual-tags-placement-1-2">function manual\_tags\_placement [1/2]</a>


```cpp
virtual void Nyx::manual_tags_placement (
    amrex::TagBoxArray & tags,
    const amrex::Vector< amrex::IntVect > & bf_lev
) override
```


Called in grid\_places after other tagging routines to modify the list of tagged points 


        

### <a href="#function-movekickdriftexact-1-2" id="function-movekickdriftexact-1-2">function moveKickDriftExact [1/2]</a>


```cpp
void Nyx::moveKickDriftExact (
    amrex::Real dt
) 
```



### <a href="#function-movekickdriftexact-1-2" id="function-movekickdriftexact-1-2">function moveKickDriftExact [1/2]</a>


```cpp
void Nyx::moveKickDriftExact (
    amrex::Real dt
) 
```



### <a href="#function-movekickexact-1-2" id="function-movekickexact-1-2">function moveKickExact [1/2]</a>


```cpp
void Nyx::moveKickExact (
    amrex::Real dt
) 
```



### <a href="#function-movekickexact-1-2" id="function-movekickexact-1-2">function moveKickExact [1/2]</a>


```cpp
void Nyx::moveKickExact (
    amrex::Real dt
) 
```



### <a href="#function-oktocontinue-1-2" id="function-oktocontinue-1-2">function okToContinue [1/2]</a>


```cpp
virtual int Nyx::okToContinue () 
```


Implements [*amrex::AmrLevel::okToContinue*](classamrex_1_1AmrLevel.md#function-oktocontinue)


### <a href="#function-oktocontinue-2-2" id="function-oktocontinue-2-2">function okToContinue [2/2]</a>


```cpp
virtual int Nyx::okToContinue () 
```


Implements [*amrex::AmrLevel::okToContinue*](classamrex_1_1AmrLevel.md#function-oktocontinue)


### <a href="#function-particle-check-point-1-2" id="function-particle-check-point-1-2">function particle\_check\_point [1/2]</a>


```cpp
void Nyx::particle_check_point (
    const std::string & dir
) 
```



### <a href="#function-particle-check-point-1-2" id="function-particle-check-point-1-2">function particle\_check\_point [1/2]</a>


```cpp
void Nyx::particle_check_point (
    const std::string & dir
) 
```



### <a href="#function-particle-derive-1-2" id="function-particle-derive-1-2">function particle\_derive [1/2]</a>


```cpp
std::unique_ptr< amrex::MultiFab > Nyx::particle_derive (
    const std::string & name,
    amrex::Real time,
    int ngrow
) 
```



### <a href="#function-particle-derive-1-2" id="function-particle-derive-1-2">function particle\_derive [1/2]</a>


```cpp
std::unique_ptr< amrex::MultiFab > Nyx::particle_derive (
    const std::string & name,
    amrex::Real time,
    int ngrow
) 
```



### <a href="#function-particle-est-time-step-1-2" id="function-particle-est-time-step-1-2">function particle\_est\_time\_step [1/2]</a>


```cpp
void Nyx::particle_est_time_step (
    amrex::Real & est_dt
) 
```



### <a href="#function-particle-est-time-step-1-2" id="function-particle-est-time-step-1-2">function particle\_est\_time\_step [1/2]</a>


```cpp
void Nyx::particle_est_time_step (
    amrex::Real & est_dt
) 
```



### <a href="#function-particle-move-random-1-2" id="function-particle-move-random-1-2">function particle\_move\_random [1/2]</a>


```cpp
void Nyx::particle_move_random () 
```



### <a href="#function-particle-move-random-1-2" id="function-particle-move-random-1-2">function particle\_move\_random [1/2]</a>


```cpp
void Nyx::particle_move_random () 
```



### <a href="#function-particle-plot-file-1-2" id="function-particle-plot-file-1-2">function particle\_plot\_file [1/2]</a>


```cpp
void Nyx::particle_plot_file (
    const std::string & dir
) 
```



### <a href="#function-particle-plot-file-1-2" id="function-particle-plot-file-1-2">function particle\_plot\_file [1/2]</a>


```cpp
void Nyx::particle_plot_file (
    const std::string & dir
) 
```



### <a href="#function-particle-post-restart-1-2" id="function-particle-post-restart-1-2">function particle\_post\_restart [1/2]</a>


```cpp
void Nyx::particle_post_restart (
    const std::string & restart_file,
    bool is_checkpoint=true
) 
```



### <a href="#function-particle-post-restart-1-2" id="function-particle-post-restart-1-2">function particle\_post\_restart [1/2]</a>


```cpp
void Nyx::particle_post_restart (
    const std::string & restart_file,
    bool is_checkpoint=true
) 
```



### <a href="#function-particle-redistribute-1-2" id="function-particle-redistribute-1-2">function particle\_redistribute [1/2]</a>


```cpp
void Nyx::particle_redistribute (
    int lbase=0,
    bool init=false
) 
```



### <a href="#function-particle-redistribute-2-2" id="function-particle-redistribute-2-2">function particle\_redistribute [2/2]</a>


```cpp
void Nyx::particle_redistribute (
    int lbase=0,
    int iteration=0,
    bool init=false
) 
```



### <a href="#function-plot-z-est-time-step-1-2" id="function-plot-z-est-time-step-1-2">function plot\_z\_est\_time\_step [1/2]</a>


```cpp
void Nyx::plot_z_est_time_step (
    amrex::Real & est_dt,
    bool & dt_changed
) 
```



### <a href="#function-plot-z-est-time-step-1-2" id="function-plot-z-est-time-step-1-2">function plot\_z\_est\_time\_step [1/2]</a>


```cpp
void Nyx::plot_z_est_time_step (
    amrex::Real & est_dt,
    bool & dt_changed
) 
```



### <a href="#function-postcoarsetimestep-1-2" id="function-postcoarsetimestep-1-2">function postCoarseTimeStep [1/2]</a>


```cpp
virtual void Nyx::postCoarseTimeStep (
    amrex::Real cumtime
) 
```



### <a href="#function-postcoarsetimestep-1-2" id="function-postcoarsetimestep-1-2">function postCoarseTimeStep [1/2]</a>


```cpp
virtual void Nyx::postCoarseTimeStep (
    amrex::Real cumtime
) 
```



### <a href="#function-post-init-1-2" id="function-post-init-1-2">function post\_init [1/2]</a>


```cpp
virtual void Nyx::post_init (
    amrex::Real stop_time
) 
```



### <a href="#function-post-init-1-2" id="function-post-init-1-2">function post\_init [1/2]</a>


```cpp
virtual void Nyx::post_init (
    amrex::Real stop_time
) 
```



### <a href="#function-post-regrid-1-2" id="function-post-regrid-1-2">function post\_regrid [1/2]</a>


```cpp
virtual void Nyx::post_regrid (
    int lbase,
    int new_finest
) 
```



### <a href="#function-post-regrid-2-2" id="function-post-regrid-2-2">function post\_regrid [2/2]</a>


```cpp
virtual void Nyx::post_regrid (
    int lbase,
    int iteration,
    int new_finest
) 
```


Implements [*amrex::AmrLevel::post\_regrid*](classamrex_1_1AmrLevel.md#function-post-regrid)


### <a href="#function-post-restart-1-2" id="function-post-restart-1-2">function post\_restart [1/2]</a>


```cpp
virtual void Nyx::post_restart () 
```


Implements [*amrex::AmrLevel::post\_restart*](classamrex_1_1AmrLevel.md#function-post-restart)


### <a href="#function-post-restart-2-2" id="function-post-restart-2-2">function post\_restart [2/2]</a>


```cpp
virtual void Nyx::post_restart () 
```


Implements [*amrex::AmrLevel::post\_restart*](classamrex_1_1AmrLevel.md#function-post-restart)


### <a href="#function-post-timestep-1-2" id="function-post-timestep-1-2">function post\_timestep [1/2]</a>


```cpp
virtual void Nyx::post_timestep (
    int iteration
) 
```


Implements [*amrex::AmrLevel::post\_timestep*](classamrex_1_1AmrLevel.md#function-post-timestep)


### <a href="#function-post-timestep-2-2" id="function-post-timestep-2-2">function post\_timestep [2/2]</a>


```cpp
virtual void Nyx::post_timestep (
    int iteration
) 
```


Implements [*amrex::AmrLevel::post\_timestep*](classamrex_1_1AmrLevel.md#function-post-timestep)


### <a href="#function-primitive-to-conserved-1-2" id="function-primitive-to-conserved-1-2">function primitive\_to\_conserved [1/2]</a>


```cpp
void Nyx::primitive_to_conserved (
    amrex::MultiFab & state
) 
```



### <a href="#function-primitive-to-conserved-1-2" id="function-primitive-to-conserved-1-2">function primitive\_to\_conserved [1/2]</a>


```cpp
void Nyx::primitive_to_conserved (
    amrex::MultiFab & state
) 
```



### <a href="#function-remove-ghost-particles-1-2" id="function-remove-ghost-particles-1-2">function remove\_ghost\_particles [1/2]</a>


```cpp
void Nyx::remove_ghost_particles () 
```



### <a href="#function-remove-ghost-particles-1-2" id="function-remove-ghost-particles-1-2">function remove\_ghost\_particles [1/2]</a>


```cpp
void Nyx::remove_ghost_particles () 
```



### <a href="#function-remove-virtual-particles-1-2" id="function-remove-virtual-particles-1-2">function remove\_virtual\_particles [1/2]</a>


```cpp
void Nyx::remove_virtual_particles () 
```



### <a href="#function-remove-virtual-particles-1-2" id="function-remove-virtual-particles-1-2">function remove\_virtual\_particles [1/2]</a>


```cpp
void Nyx::remove_virtual_particles () 
```



### <a href="#function-reset-internal-energy-1-2" id="function-reset-internal-energy-1-2">function reset\_internal\_energy [1/2]</a>


```cpp
void Nyx::reset_internal_energy (
    amrex::MultiFab & State,
    amrex::MultiFab & DiagEOS,
    amrex::MultiFab & reset_e_src
) 
```



### <a href="#function-reset-internal-energy-1-2" id="function-reset-internal-energy-1-2">function reset\_internal\_energy [1/2]</a>


```cpp
void Nyx::reset_internal_energy (
    amrex::MultiFab & State,
    amrex::MultiFab & DiagEOS,
    amrex::MultiFab & reset_e_src
) 
```



### <a href="#function-restart-1-2" id="function-restart-1-2">function restart [1/2]</a>


```cpp
virtual void Nyx::restart (
    amrex::Amr & papa,
    istream & is,
    bool b_read_special=false
) 
```



### <a href="#function-restart-2-2" id="function-restart-2-2">function restart [2/2]</a>


```cpp
virtual void Nyx::restart (
    amrex::Amr & papa,
    istream & is,
    bool b_read_special=false
) 
```



### <a href="#function-setplotvariables-1-2" id="function-setplotvariables-1-2">function setPlotVariables [1/2]</a>


```cpp
virtual void Nyx::setPlotVariables () 
```


Implements [*amrex::AmrLevel::setPlotVariables*](classamrex_1_1AmrLevel.md#function-setplotvariables)


### <a href="#function-setplotvariables-2-2" id="function-setplotvariables-2-2">function setPlotVariables [2/2]</a>


```cpp
virtual void Nyx::setPlotVariables () 
```


Implements [*amrex::AmrLevel::setPlotVariables*](classamrex_1_1AmrLevel.md#function-setplotvariables)


### <a href="#function-settimelevel-1-2" id="function-settimelevel-1-2">function setTimeLevel [1/2]</a>


```cpp
virtual void Nyx::setTimeLevel (
    amrex::Real time,
    amrex::Real dt_old,
    amrex::Real dt_new
) 
```



### <a href="#function-settimelevel-1-2" id="function-settimelevel-1-2">function setTimeLevel [1/2]</a>


```cpp
virtual void Nyx::setTimeLevel (
    amrex::Real time,
    amrex::Real dt_old,
    amrex::Real dt_new
) 
```



### <a href="#function-setup-ghost-particles-1-2" id="function-setup-ghost-particles-1-2">function setup\_ghost\_particles [1/2]</a>


```cpp
void Nyx::setup_ghost_particles (
    int ngrow
) 
```



### <a href="#function-setup-ghost-particles-1-2" id="function-setup-ghost-particles-1-2">function setup\_ghost\_particles [1/2]</a>


```cpp
void Nyx::setup_ghost_particles (
    int ngrow
) 
```



### <a href="#function-setup-virtual-particles-1-2" id="function-setup-virtual-particles-1-2">function setup\_virtual\_particles [1/2]</a>


```cpp
void Nyx::setup_virtual_particles () 
```



### <a href="#function-setup-virtual-particles-1-2" id="function-setup-virtual-particles-1-2">function setup\_virtual\_particles [1/2]</a>


```cpp
void Nyx::setup_virtual_particles () 
```



### <a href="#function-strang-first-step-1-2" id="function-strang-first-step-1-2">function strang\_first\_step [1/2]</a>


```cpp
void Nyx::strang_first_step (
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab & state,
    amrex::MultiFab & dstate
) 
```



### <a href="#function-strang-first-step-1-2" id="function-strang-first-step-1-2">function strang\_first\_step [1/2]</a>


```cpp
void Nyx::strang_first_step (
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab & state,
    amrex::MultiFab & dstate
) 
```



### <a href="#function-strang-hydro-1-2" id="function-strang-hydro-1-2">function strang\_hydro [1/2]</a>


```cpp
void Nyx::strang_hydro (
    amrex::Real time,
    amrex::Real dt,
    amrex::Real a_old,
    amrex::Real a_new
) 
```



### <a href="#function-strang-hydro-1-2" id="function-strang-hydro-1-2">function strang\_hydro [1/2]</a>


```cpp
void Nyx::strang_hydro (
    amrex::Real time,
    amrex::Real dt,
    amrex::Real a_old,
    amrex::Real a_new
) 
```



### <a href="#function-strang-second-step-1-2" id="function-strang-second-step-1-2">function strang\_second\_step [1/2]</a>


```cpp
void Nyx::strang_second_step (
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab & state,
    amrex::MultiFab & dstate
) 
```



### <a href="#function-strang-second-step-1-2" id="function-strang-second-step-1-2">function strang\_second\_step [1/2]</a>


```cpp
void Nyx::strang_second_step (
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab & state,
    amrex::MultiFab & dstate
) 
```



### <a href="#function-theplotfiletype-1-2" id="function-theplotfiletype-1-2">function thePlotFileType [1/2]</a>


```cpp
virtual std::string Nyx::thePlotFileType () const
```


A string written as the first item in `write_plot_file()` at level zero. It is so we can distinguish between different types of plot files. For [**Nyx**](classNyx.md) it has the form: Nyx-Vnnn. 


        
Implements [*amrex::AmrLevel::thePlotFileType*](classamrex_1_1AmrLevel.md#function-theplotfiletype)


### <a href="#function-theplotfiletype-2-2" id="function-theplotfiletype-2-2">function thePlotFileType [2/2]</a>


```cpp
virtual std::string Nyx::thePlotFileType () const
```


A string written as the first item in `write_plot_file()` at level zero. It is so we can distinguish between different types of plot files. For [**Nyx**](classNyx.md) it has the form: Nyx-Vnnn. 


        
Implements [*amrex::AmrLevel::thePlotFileType*](classamrex_1_1AmrLevel.md#function-theplotfiletype)


### <a href="#function-time-center-source-terms-1-2" id="function-time-center-source-terms-1-2">function time\_center\_source\_terms [1/2]</a>


```cpp
void Nyx::time_center_source_terms (
    amrex::MultiFab & S_new,
    amrex::MultiFab & ext_src_old,
    amrex::MultiFab & ext_src_new,
    amrex::Real dt
) 
```



### <a href="#function-time-center-source-terms-1-2" id="function-time-center-source-terms-1-2">function time\_center\_source\_terms [1/2]</a>


```cpp
void Nyx::time_center_source_terms (
    amrex::MultiFab & S_new,
    amrex::MultiFab & ext_src_old,
    amrex::MultiFab & ext_src_new,
    amrex::Real dt
) 
```



### <a href="#function-update-state-with-sources-1-2" id="function-update-state-with-sources-1-2">function update\_state\_with\_sources [1/2]</a>


```cpp
void Nyx::update_state_with_sources (
    amrex::MultiFab & S_old,
    amrex::MultiFab & S_new,
    amrex::MultiFab & ext_src_old,
    amrex::MultiFab & hydro_src,
    amrex::MultiFab & grav,
    amrex::MultiFab & divu_cc,
    amrex::Real dt,
    amrex::Real a_old,
    amrex::Real a_new
) 
```



### <a href="#function-update-state-with-sources-1-2" id="function-update-state-with-sources-1-2">function update\_state\_with\_sources [1/2]</a>


```cpp
void Nyx::update_state_with_sources (
    amrex::MultiFab & S_old,
    amrex::MultiFab & S_new,
    amrex::MultiFab & ext_src_old,
    amrex::MultiFab & hydro_src,
    amrex::MultiFab & grav,
    amrex::MultiFab & divu_cc,
    amrex::Real dt,
    amrex::Real a_old,
    amrex::Real a_new
) 
```



### <a href="#function-vol-weight-squared-sum-1-2" id="function-vol-weight-squared-sum-1-2">function vol\_weight\_squared\_sum [1/2]</a>


```cpp
amrex::Real Nyx::vol_weight_squared_sum (
    const std::string & name,
    amrex::Real time
) 
```



### <a href="#function-vol-weight-squared-sum-1-2" id="function-vol-weight-squared-sum-1-2">function vol\_weight\_squared\_sum [1/2]</a>


```cpp
amrex::Real Nyx::vol_weight_squared_sum (
    const std::string & name,
    amrex::Real time
) 
```



### <a href="#function-vol-weight-squared-sum-level-1-2" id="function-vol-weight-squared-sum-level-1-2">function vol\_weight\_squared\_sum\_level [1/2]</a>


```cpp
amrex::Real Nyx::vol_weight_squared_sum_level (
    const std::string & name,
    amrex::Real time
) 
```



### <a href="#function-vol-weight-squared-sum-level-1-2" id="function-vol-weight-squared-sum-level-1-2">function vol\_weight\_squared\_sum\_level [1/2]</a>


```cpp
amrex::Real Nyx::vol_weight_squared_sum_level (
    const std::string & name,
    amrex::Real time
) 
```



### <a href="#function-vol-weight-sum-1-4" id="function-vol-weight-sum-1-4">function vol\_weight\_sum [1/4]</a>


```cpp
amrex::Real Nyx::vol_weight_sum (
    const std::string & name,
    amrex::Real time,
    bool masked
) 
```



### <a href="#function-vol-weight-sum-2-4" id="function-vol-weight-sum-2-4">function vol\_weight\_sum [2/4]</a>


```cpp
amrex::Real Nyx::vol_weight_sum (
    amrex::MultiFab & mf,
    bool masked
) 
```



### <a href="#function-vol-weight-sum-1-4" id="function-vol-weight-sum-1-4">function vol\_weight\_sum [1/4]</a>


```cpp
amrex::Real Nyx::vol_weight_sum (
    const std::string & name,
    amrex::Real time,
    bool masked
) 
```



### <a href="#function-vol-weight-sum-2-4" id="function-vol-weight-sum-2-4">function vol\_weight\_sum [2/4]</a>


```cpp
amrex::Real Nyx::vol_weight_sum (
    amrex::MultiFab & mf,
    bool masked
) 
```



### <a href="#function-writejobinfo-1-2" id="function-writejobinfo-1-2">function writeJobInfo [1/2]</a>


```cpp
void Nyx::writeJobInfo (
    const std::string & dir
) 
```



### <a href="#function-writejobinfo-1-2" id="function-writejobinfo-1-2">function writeJobInfo [1/2]</a>


```cpp
void Nyx::writeJobInfo (
    const std::string & dir
) 
```



### <a href="#function-writemultifabasplotfile-1-2" id="function-writemultifabasplotfile-1-2">function writeMultiFabAsPlotFile [1/2]</a>


```cpp
void Nyx::writeMultiFabAsPlotFile (
    const std::string & pltfile,
    const amrex::MultiFab & mf,
    std::string componentName
) 
```



### <a href="#function-writemultifabasplotfile-1-2" id="function-writemultifabasplotfile-1-2">function writeMultiFabAsPlotFile [1/2]</a>


```cpp
void Nyx::writeMultiFabAsPlotFile (
    const std::string & pltfile,
    const amrex::MultiFab & mf,
    std::string componentName
) 
```



### <a href="#function-writeplotfile-1-2" id="function-writeplotfile-1-2">function writePlotFile [1/2]</a>


```cpp
virtual void Nyx::writePlotFile (
    const std::string & dir,
    ostream & os,
    amrex::VisMF::How how
) 
```



### <a href="#function-writeplotfile-1-2" id="function-writeplotfile-1-2">function writePlotFile [1/2]</a>


```cpp
virtual void Nyx::writePlotFile (
    const std::string & dir,
    ostream & os,
    amrex::VisMF::How how
) 
```



### <a href="#function-writeplotfilepost-1-2" id="function-writeplotfilepost-1-2">function writePlotFilePost [1/2]</a>


```cpp
virtual void Nyx::writePlotFilePost (
    const std::string & dir,
    ostream & os
) 
```



### <a href="#function-writeplotfilepost-2-2" id="function-writeplotfilepost-2-2">function writePlotFilePost [2/2]</a>


```cpp
virtual void Nyx::writePlotFilePost (
    const std::string & dir,
    ostream & os
) 
```



### <a href="#function-writeplotfilepre-1-2" id="function-writeplotfilepre-1-2">function writePlotFilePre [1/2]</a>


```cpp
virtual void Nyx::writePlotFilePre (
    const std::string & dir,
    ostream & os
) 
```



### <a href="#function-writeplotfilepre-2-2" id="function-writeplotfilepre-2-2">function writePlotFilePre [2/2]</a>


```cpp
virtual void Nyx::writePlotFilePre (
    const std::string & dir,
    ostream & os
) 
```



### <a href="#function-writeplotnow-1-2" id="function-writeplotnow-1-2">function writePlotNow [1/2]</a>


```cpp
virtual bool Nyx::writePlotNow () 
```


Implements [*amrex::AmrLevel::writePlotNow*](classamrex_1_1AmrLevel.md#function-writeplotnow)


### <a href="#function-writeplotnow-1-2" id="function-writeplotnow-1-2">function writePlotNow [1/2]</a>


```cpp
virtual bool Nyx::writePlotNow () 
```


Implements [*amrex::AmrLevel::writePlotNow*](classamrex_1_1AmrLevel.md#function-writeplotnow)


### <a href="#function-write-parameter-file-1-2" id="function-write-parameter-file-1-2">function write\_parameter\_file [1/2]</a>


```cpp
virtual void Nyx::write_parameter_file (
    const std::string & dir
) 
```



### <a href="#function-write-parameter-file-2-2" id="function-write-parameter-file-2-2">function write\_parameter\_file [2/2]</a>


```cpp
virtual void Nyx::write_parameter_file (
    const std::string & dir
) 
```



### <a href="#function-nyx-1-2" id="function-nyx-1-2">function ~Nyx [1/2]</a>


```cpp
virtual Nyx::~Nyx () 
```



### <a href="#function-nyx-2-2" id="function-nyx-2-2">function ~Nyx [2/2]</a>


```cpp
virtual Nyx::~Nyx () 
```


## Public Static Functions Documentation


### <a href="#function-do-hydro-1-2" id="function-do-hydro-1-2">function Do\_Hydro [1/2]</a>


```cpp
static inline int Nyx::Do_Hydro () 
```



### <a href="#function-do-hydro-2-2" id="function-do-hydro-2-2">function Do\_Hydro [2/2]</a>


```cpp
static int Nyx::Do_Hydro () 
```



### <a href="#function-initderivelist-1-2" id="function-initderivelist-1-2">function InitDeriveList [1/2]</a>


```cpp
static void Nyx::InitDeriveList () 
```



### <a href="#function-initderivelist-2-2" id="function-initderivelist-2-2">function InitDeriveList [2/2]</a>


```cpp
static void Nyx::InitDeriveList () 
```



### <a href="#function-initerrorlist-1-2" id="function-initerrorlist-1-2">function InitErrorList [1/2]</a>


```cpp
static void Nyx::InitErrorList () 
```



### <a href="#function-initerrorlist-2-2" id="function-initerrorlist-2-2">function InitErrorList [2/2]</a>


```cpp
static void Nyx::InitErrorList () 
```



### <a href="#function-alloc-simd-vec-1-2" id="function-alloc-simd-vec-1-2">function alloc\_simd\_vec [1/2]</a>


```cpp
static void Nyx::alloc_simd_vec () 
```



### <a href="#function-alloc-simd-vec-1-2" id="function-alloc-simd-vec-1-2">function alloc\_simd\_vec [1/2]</a>


```cpp
static void Nyx::alloc_simd_vec () 
```



### <a href="#function-dealloc-simd-vec-1-2" id="function-dealloc-simd-vec-1-2">function dealloc\_simd\_vec [1/2]</a>


```cpp
static void Nyx::dealloc_simd_vec () 
```



### <a href="#function-dealloc-simd-vec-1-2" id="function-dealloc-simd-vec-1-2">function dealloc\_simd\_vec [1/2]</a>


```cpp
static void Nyx::dealloc_simd_vec () 
```



### <a href="#function-error-setup-1-2" id="function-error-setup-1-2">function error\_setup [1/2]</a>


```cpp
static void Nyx::error_setup () 
```



### <a href="#function-error-setup-2-2" id="function-error-setup-2-2">function error\_setup [2/2]</a>


```cpp
static void Nyx::error_setup () 
```



### <a href="#function-hydro-setup-1-2" id="function-hydro-setup-1-2">function hydro\_setup [1/2]</a>


```cpp
static void Nyx::hydro_setup () 
```



### <a href="#function-hydro-setup-2-2" id="function-hydro-setup-2-2">function hydro\_setup [2/2]</a>


```cpp
static void Nyx::hydro_setup () 
```



### <a href="#function-no-hydro-setup-1-2" id="function-no-hydro-setup-1-2">function no\_hydro\_setup [1/2]</a>


```cpp
static void Nyx::no_hydro_setup () 
```



### <a href="#function-no-hydro-setup-1-2" id="function-no-hydro-setup-1-2">function no\_hydro\_setup [1/2]</a>


```cpp
static void Nyx::no_hydro_setup () 
```



### <a href="#function-num-grow-1-2" id="function-num-grow-1-2">function num\_grow [1/2]</a>


```cpp
static inline int Nyx::num_grow () 
```



### <a href="#function-num-grow-2-2" id="function-num-grow-2-2">function num\_grow [2/2]</a>


```cpp
static int Nyx::num_grow () 
```



### <a href="#function-read-comoving-params-1-2" id="function-read-comoving-params-1-2">function read\_comoving\_params [1/2]</a>


```cpp
static void Nyx::read_comoving_params () 
```



### <a href="#function-read-comoving-params-2-2" id="function-read-comoving-params-2-2">function read\_comoving\_params [2/2]</a>


```cpp
static void Nyx::read_comoving_params () 
```



### <a href="#function-read-init-params-1-2" id="function-read-init-params-1-2">function read\_init\_params [1/2]</a>


```cpp
static void Nyx::read_init_params () 
```



### <a href="#function-read-init-params-2-2" id="function-read-init-params-2-2">function read\_init\_params [2/2]</a>


```cpp
static void Nyx::read_init_params () 
```



### <a href="#function-read-particle-params-1-2" id="function-read-particle-params-1-2">function read\_particle\_params [1/2]</a>


```cpp
static void Nyx::read_particle_params () 
```



### <a href="#function-read-particle-params-2-2" id="function-read-particle-params-2-2">function read\_particle\_params [2/2]</a>


```cpp
static void Nyx::read_particle_params () 
```



### <a href="#function-set-simd-width-1-2" id="function-set-simd-width-1-2">function set\_simd\_width [1/2]</a>


```cpp
static void Nyx::set_simd_width (
    const int simd_width
) 
```



### <a href="#function-set-simd-width-1-2" id="function-set-simd-width-1-2">function set\_simd\_width [1/2]</a>


```cpp
static void Nyx::set_simd_width (
    const int simd_width
) 
```



### <a href="#function-theactiveparticles-1-2" id="function-theactiveparticles-1-2">function theActiveParticles [1/2]</a>


```cpp
static amrex::Vector< NyxParticleContainerBase * > & Nyx::theActiveParticles () 
```



### <a href="#function-theactiveparticles-2-2" id="function-theactiveparticles-2-2">function theActiveParticles [2/2]</a>


```cpp
static amrex::Vector< NyxParticleContainerBase * > & Nyx::theActiveParticles () 
```



### <a href="#function-thedmpc-1-2" id="function-thedmpc-1-2">function theDMPC [1/2]</a>


```cpp
static DarkMatterParticleContainer * Nyx::theDMPC () 
```



### <a href="#function-thedmpc-2-2" id="function-thedmpc-2-2">function theDMPC [2/2]</a>


```cpp
static DarkMatterParticleContainer * Nyx::theDMPC () 
```



### <a href="#function-theghostpc-1-2" id="function-theghostpc-1-2">function theGhostPC [1/2]</a>


```cpp
static DarkMatterParticleContainer * Nyx::theGhostPC () 
```



### <a href="#function-theghostpc-2-2" id="function-theghostpc-2-2">function theGhostPC [2/2]</a>


```cpp
static DarkMatterParticleContainer * Nyx::theGhostPC () 
```



### <a href="#function-theghostparticles-1-2" id="function-theghostparticles-1-2">function theGhostParticles [1/2]</a>


```cpp
static amrex::Vector< NyxParticleContainerBase * > & Nyx::theGhostParticles () 
```



### <a href="#function-theghostparticles-2-2" id="function-theghostparticles-2-2">function theGhostParticles [2/2]</a>


```cpp
static amrex::Vector< NyxParticleContainerBase * > & Nyx::theGhostParticles () 
```



### <a href="#function-theghostspc-1-2" id="function-theghostspc-1-2">function theGhostSPC [1/2]</a>


```cpp
static StellarParticleContainer * Nyx::theGhostSPC () 
```



### <a href="#function-theghostspc-2-2" id="function-theghostspc-2-2">function theGhostSPC [2/2]</a>


```cpp
static StellarParticleContainer * Nyx::theGhostSPC () 
```



### <a href="#function-thespc-1-2" id="function-thespc-1-2">function theSPC [1/2]</a>


```cpp
static StellarParticleContainer * Nyx::theSPC () 
```



### <a href="#function-thespc-2-2" id="function-thespc-2-2">function theSPC [2/2]</a>


```cpp
static StellarParticleContainer * Nyx::theSPC () 
```



### <a href="#function-thevirtpc-1-2" id="function-thevirtpc-1-2">function theVirtPC [1/2]</a>


```cpp
static DarkMatterParticleContainer * Nyx::theVirtPC () 
```



### <a href="#function-thevirtpc-2-2" id="function-thevirtpc-2-2">function theVirtPC [2/2]</a>


```cpp
static DarkMatterParticleContainer * Nyx::theVirtPC () 
```



### <a href="#function-thevirtspc-1-2" id="function-thevirtspc-1-2">function theVirtSPC [1/2]</a>


```cpp
static StellarParticleContainer * Nyx::theVirtSPC () 
```



### <a href="#function-thevirtspc-2-2" id="function-thevirtspc-2-2">function theVirtSPC [2/2]</a>


```cpp
static StellarParticleContainer * Nyx::theVirtSPC () 
```



### <a href="#function-thevirtualparticles-1-2" id="function-thevirtualparticles-1-2">function theVirtualParticles [1/2]</a>


```cpp
static amrex::Vector< NyxParticleContainerBase * > & Nyx::theVirtualParticles () 
```



### <a href="#function-thevirtualparticles-2-2" id="function-thevirtualparticles-2-2">function theVirtualParticles [2/2]</a>


```cpp
static amrex::Vector< NyxParticleContainerBase * > & Nyx::theVirtualParticles () 
```



### <a href="#function-variable-cleanup-1-2" id="function-variable-cleanup-1-2">function variable\_cleanup [1/2]</a>


```cpp
static void Nyx::variable_cleanup () 
```



### <a href="#function-variable-cleanup-2-2" id="function-variable-cleanup-2-2">function variable\_cleanup [2/2]</a>


```cpp
static void Nyx::variable_cleanup () 
```



### <a href="#function-variable-setup-1-2" id="function-variable-setup-1-2">function variable\_setup [1/2]</a>


```cpp
static void Nyx::variable_setup () 
```



### <a href="#function-variable-setup-2-2" id="function-variable-setup-2-2">function variable\_setup [2/2]</a>


```cpp
static void Nyx::variable_setup () 
```


## Protected Attributes Documentation


### <a href="#variable-fillpatchedoldstate-ok" id="variable-fillpatchedoldstate-ok">variable FillPatchedOldState\_ok </a>


```cpp
bool Nyx::FillPatchedOldState_ok;
```



### <a href="#variable-flux-reg" id="variable-flux-reg">variable flux\_reg </a>


```cpp
amrex::FluxRegister * Nyx::flux_reg;
```


## Protected Static Attributes Documentation


### <a href="#variable-num-grow" id="variable-num-grow">variable NUM\_GROW </a>


```cpp
static int Nyx::NUM_GROW;
```



### <a href="#variable-nrep" id="variable-nrep">variable Nrep </a>


```cpp
IntVect Nyx::Nrep;
```



### <a href="#variable-add-ext-src" id="variable-add-ext-src">variable add\_ext\_src </a>


```cpp
static int Nyx::add_ext_src;
```



### <a href="#variable-allow-untagging" id="variable-allow-untagging">variable allow\_untagging </a>


```cpp
static int Nyx::allow_untagging;
```



### <a href="#variable-analysis-z-values" id="variable-analysis-z-values">variable analysis\_z\_values </a>


```cpp
static amrex::Vector< amrex::Real > Nyx::analysis_z_values;
```



### <a href="#variable-average-dm-density" id="variable-average-dm-density">variable average\_dm\_density </a>


```cpp
static amrex::Real Nyx::average_dm_density;
```



### <a href="#variable-average-gas-density" id="variable-average-gas-density">variable average\_gas\_density </a>


```cpp
static amrex::Real Nyx::average_gas_density;
```



### <a href="#variable-average-neutr-density" id="variable-average-neutr-density">variable average\_neutr\_density </a>


```cpp
static amrex::Real Nyx::average_neutr_density;
```



### <a href="#variable-average-total-density" id="variable-average-total-density">variable average\_total\_density </a>


```cpp
static amrex::Real Nyx::average_total_density;
```



### <a href="#variable-cfl" id="variable-cfl">variable cfl </a>


```cpp
static amrex::Real Nyx::cfl;
```



### <a href="#variable-change-max" id="variable-change-max">variable change\_max </a>


```cpp
static amrex::Real Nyx::change_max;
```



### <a href="#variable-corner-coupling" id="variable-corner-coupling">variable corner\_coupling </a>


```cpp
static int Nyx::corner_coupling;
```



### <a href="#variable-do-dm-particles" id="variable-do-dm-particles">variable do\_dm\_particles </a>


```cpp
bool Nyx::do_dm_particles;
```



### <a href="#variable-do-forcing" id="variable-do-forcing">variable do\_forcing </a>


```cpp
static int Nyx::do_forcing;
```



### <a href="#variable-do-grav" id="variable-do-grav">variable do\_grav </a>


```cpp
static int Nyx::do_grav;
```



### <a href="#variable-do-hydro" id="variable-do-hydro">variable do\_hydro </a>


```cpp
static int Nyx::do_hydro;
```



### <a href="#variable-do-reflux" id="variable-do-reflux">variable do\_reflux </a>


```cpp
static int Nyx::do_reflux;
```



### <a href="#variable-do-special-tagging" id="variable-do-special-tagging">variable do\_special\_tagging </a>


```cpp
static int Nyx::do_special_tagging;
```



### <a href="#variable-dump-old" id="variable-dump-old">variable dump\_old </a>


```cpp
static bool Nyx::dump_old;
```



### <a href="#variable-err-list" id="variable-err-list">variable err\_list </a>


```cpp
static amrex::ErrorList Nyx::err_list;
```



### <a href="#variable-gamma" id="variable-gamma">variable gamma </a>


```cpp
static amrex::Real Nyx::gamma;
```



### <a href="#variable-h-species" id="variable-h-species">variable h\_species </a>


```cpp
static amrex::Real Nyx::h_species;
```



### <a href="#variable-he-species" id="variable-he-species">variable he\_species </a>


```cpp
static amrex::Real Nyx::he_species;
```



### <a href="#variable-heat-cool-type" id="variable-heat-cool-type">variable heat\_cool\_type </a>


```cpp
static int Nyx::heat_cool_type;
```



### <a href="#variable-inhomo-grid" id="variable-inhomo-grid">variable inhomo\_grid </a>


```cpp
static int Nyx::inhomo_grid;
```



### <a href="#variable-inhomo-reion" id="variable-inhomo-reion">variable inhomo\_reion </a>


```cpp
static int Nyx::inhomo_reion;
```



### <a href="#variable-inhomo-zhi-file" id="variable-inhomo-zhi-file">variable inhomo\_zhi\_file </a>


```cpp
static std::string Nyx::inhomo_zhi_file;
```



### <a href="#variable-init-shrink" id="variable-init-shrink">variable init\_shrink </a>


```cpp
static amrex::Real Nyx::init_shrink;
```



### <a href="#variable-normalize-species" id="variable-normalize-species">variable normalize\_species </a>


```cpp
static int Nyx::normalize_species;
```



### <a href="#variable-nsteps-from-plotfile" id="variable-nsteps-from-plotfile">variable nsteps\_from\_plotfile </a>


```cpp
static int Nyx::nsteps_from_plotfile;
```



### <a href="#variable-num-particle-ghosts" id="variable-num-particle-ghosts">variable num\_particle\_ghosts </a>


```cpp
int Nyx::num_particle_ghosts;
```



### <a href="#variable-particle-init-type" id="variable-particle-init-type">variable particle\_init\_type </a>


```cpp
std::string Nyx::particle_init_type;
```


How do we want to initialize the particles? Must be "Random", "Cosmological" or "AsciiFile" 


        

### <a href="#variable-particle-initrandom-count" id="variable-particle-initrandom-count">variable particle\_initrandom\_count </a>


```cpp
long Nyx::particle_initrandom_count;
```



### <a href="#variable-particle-initrandom-count-per-box" id="variable-particle-initrandom-count-per-box">variable particle\_initrandom\_count\_per\_box </a>


```cpp
long Nyx::particle_initrandom_count_per_box;
```



### <a href="#variable-particle-initrandom-iseed" id="variable-particle-initrandom-iseed">variable particle\_initrandom\_iseed </a>


```cpp
int Nyx::particle_initrandom_iseed;
```



### <a href="#variable-particle-initrandom-mass" id="variable-particle-initrandom-mass">variable particle\_initrandom\_mass </a>


```cpp
Real Nyx::particle_initrandom_mass;
```



### <a href="#variable-particle-initrandom-serialize" id="variable-particle-initrandom-serialize">variable particle\_initrandom\_serialize </a>


```cpp
bool Nyx::particle_initrandom_serialize;
```



### <a href="#variable-particle-move-type" id="variable-particle-move-type">variable particle\_move\_type </a>


```cpp
std::string Nyx::particle_move_type;
```


How do we want to move the particles? Must be "Random" or "Gravitational" 


        

### <a href="#variable-particle-skip-factor" id="variable-particle-skip-factor">variable particle\_skip\_factor </a>


```cpp
int Nyx::particle_skip_factor;
```



### <a href="#variable-phys-bc" id="variable-phys-bc">variable phys\_bc </a>


```cpp
static amrex::BCRec Nyx::phys_bc;
```



### <a href="#variable-plot-z-values" id="variable-plot-z-values">variable plot\_z\_values </a>


```cpp
static amrex::Vector< amrex::Real > Nyx::plot_z_values;
```



### <a href="#variable-ppm-flatten-before-integrals" id="variable-ppm-flatten-before-integrals">variable ppm\_flatten\_before\_integrals </a>


```cpp
static int Nyx::ppm_flatten_before_integrals;
```



### <a href="#variable-ppm-reference" id="variable-ppm-reference">variable ppm\_reference </a>


```cpp
static int Nyx::ppm_reference;
```



### <a href="#variable-ppm-type" id="variable-ppm-type">variable ppm\_type </a>


```cpp
static int Nyx::ppm_type;
```



### <a href="#variable-previouscputimeused" id="variable-previouscputimeused">variable previousCPUTimeUsed </a>


```cpp
static amrex::Real Nyx::previousCPUTimeUsed;
```


for keeping track of the amount of CPU time used  this will persist after restarts 


        

### <a href="#variable-small-dens" id="variable-small-dens">variable small\_dens </a>


```cpp
static amrex::Real Nyx::small_dens;
```



### <a href="#variable-small-temp" id="variable-small-temp">variable small\_temp </a>


```cpp
static amrex::Real Nyx::small_temp;
```



### <a href="#variable-startcputime" id="variable-startcputime">variable startCPUTime </a>


```cpp
static amrex::Real Nyx::startCPUTime;
```



### <a href="#variable-strang-split" id="variable-strang-split">variable strang\_split </a>


```cpp
static int Nyx::strang_split;
```



### <a href="#variable-use-colglaz" id="variable-use-colglaz">variable use\_colglaz </a>


```cpp
static int Nyx::use_colglaz;
```



### <a href="#variable-use-const-species" id="variable-use-const-species">variable use\_const\_species </a>


```cpp
static int Nyx::use_const_species;
```



### <a href="#variable-use-exact-gravity" id="variable-use-exact-gravity">variable use\_exact\_gravity </a>


```cpp
static int Nyx::use_exact_gravity;
```



### <a href="#variable-use-flattening" id="variable-use-flattening">variable use\_flattening </a>


```cpp
static int Nyx::use_flattening;
```



### <a href="#variable-verbose" id="variable-verbose">variable verbose </a>


```cpp
static int Nyx::verbose;
```



### <a href="#variable-version-2" id="variable-version-2">variable version\_2 </a>


```cpp
static int Nyx::version_2;
```


## Protected Functions Documentation


### <a href="#function-average-down-1-4" id="function-average-down-1-4">function average\_down [1/4]</a>


```cpp
void Nyx::average_down () 
```



### <a href="#function-average-down-2-4" id="function-average-down-2-4">function average\_down [2/4]</a>


```cpp
void Nyx::average_down (
    int state_indx
) 
```



### <a href="#function-average-down-1-4" id="function-average-down-1-4">function average\_down [1/4]</a>


```cpp
void Nyx::average_down () 
```



### <a href="#function-average-down-2-4" id="function-average-down-2-4">function average\_down [2/4]</a>


```cpp
void Nyx::average_down (
    int state_indx
) 
```



### <a href="#function-build-metrics-1-2" id="function-build-metrics-1-2">function build\_metrics [1/2]</a>


```cpp
void Nyx::build_metrics () 
```



### <a href="#function-build-metrics-1-2" id="function-build-metrics-1-2">function build\_metrics [1/2]</a>


```cpp
void Nyx::build_metrics () 
```



### <a href="#function-compute-average-density-1-2" id="function-compute-average-density-1-2">function compute\_average\_density [1/2]</a>


```cpp
void Nyx::compute_average_density () 
```



### <a href="#function-compute-average-density-1-2" id="function-compute-average-density-1-2">function compute\_average\_density [1/2]</a>


```cpp
void Nyx::compute_average_density () 
```



### <a href="#function-compute-average-species-1-2" id="function-compute-average-species-1-2">function compute\_average\_species [1/2]</a>


```cpp
void Nyx::compute_average_species (
    int nspec,
    int naux,
    amrex::Vector< amrex::Real > & average_species
) 
```



### <a href="#function-compute-average-species-1-2" id="function-compute-average-species-1-2">function compute\_average\_species [1/2]</a>


```cpp
void Nyx::compute_average_species (
    int nspec,
    int naux,
    amrex::Vector< amrex::Real > & average_species
) 
```



### <a href="#function-compute-average-temperature-1-2" id="function-compute-average-temperature-1-2">function compute\_average\_temperature [1/2]</a>


```cpp
void Nyx::compute_average_temperature (
    amrex::Real & average_temperature
) 
```



### <a href="#function-compute-average-temperature-1-2" id="function-compute-average-temperature-1-2">function compute\_average\_temperature [1/2]</a>


```cpp
void Nyx::compute_average_temperature (
    amrex::Real & average_temperature
) 
```



### <a href="#function-enforce-consistent-e-1-2" id="function-enforce-consistent-e-1-2">function enforce\_consistent\_e [1/2]</a>


```cpp
void Nyx::enforce_consistent_e (
    amrex::MultiFab & S
) 
```



### <a href="#function-enforce-consistent-e-1-2" id="function-enforce-consistent-e-1-2">function enforce\_consistent\_e [1/2]</a>


```cpp
void Nyx::enforce_consistent_e (
    amrex::MultiFab & S
) 
```



### <a href="#function-enforce-nonnegative-species-1-2" id="function-enforce-nonnegative-species-1-2">function enforce\_nonnegative\_species [1/2]</a>


```cpp
void Nyx::enforce_nonnegative_species (
    amrex::MultiFab & S_new
) 
```



### <a href="#function-enforce-nonnegative-species-1-2" id="function-enforce-nonnegative-species-1-2">function enforce\_nonnegative\_species [1/2]</a>


```cpp
void Nyx::enforce_nonnegative_species (
    amrex::MultiFab & S_new
) 
```



### <a href="#function-get-flux-reg-1-4" id="function-get-flux-reg-1-4">function get\_flux\_reg [1/4]</a>


```cpp
inline amrex::FluxRegister & Nyx::get_flux_reg () 
```



### <a href="#function-get-flux-reg-2-4" id="function-get-flux-reg-2-4">function get\_flux\_reg [2/4]</a>


```cpp
inline amrex::FluxRegister & Nyx::get_flux_reg (
    int lev
) 
```



### <a href="#function-get-flux-reg-3-4" id="function-get-flux-reg-3-4">function get\_flux\_reg [3/4]</a>


```cpp
amrex::FluxRegister & Nyx::get_flux_reg () 
```



### <a href="#function-get-flux-reg-4-4" id="function-get-flux-reg-4-4">function get\_flux\_reg [4/4]</a>


```cpp
amrex::FluxRegister & Nyx::get_flux_reg (
    int lev
) 
```



### <a href="#function-get-level-1-2" id="function-get-level-1-2">function get\_level [1/2]</a>


```cpp
inline Nyx & Nyx::get_level (
    int lev
) 
```



### <a href="#function-get-level-2-2" id="function-get-level-2-2">function get\_level [2/2]</a>


```cpp
Nyx & Nyx::get_level (
    int lev
) 
```



### <a href="#function-reflux-1-2" id="function-reflux-1-2">function reflux [1/2]</a>


```cpp
void Nyx::reflux () 
```



### <a href="#function-reflux-1-2" id="function-reflux-1-2">function reflux [1/2]</a>


```cpp
void Nyx::reflux () 
```



### <a href="#function-retrievedm-1-2" id="function-retrievedm-1-2">function retrieveDM [1/2]</a>


```cpp
std::string Nyx::retrieveDM () 
```



### <a href="#function-retrievedm-1-2" id="function-retrievedm-1-2">function retrieveDM [1/2]</a>


```cpp
std::string Nyx::retrieveDM () 
```



### <a href="#function-set-small-values-1-2" id="function-set-small-values-1-2">function set\_small\_values [1/2]</a>


```cpp
void Nyx::set_small_values () 
```



### <a href="#function-set-small-values-1-2" id="function-set-small-values-1-2">function set\_small\_values [1/2]</a>


```cpp
void Nyx::set_small_values () 
```



### <a href="#function-sum-integrated-quantities-1-2" id="function-sum-integrated-quantities-1-2">function sum\_integrated\_quantities [1/2]</a>


```cpp
virtual void Nyx::sum_integrated_quantities () 
```



### <a href="#function-sum-integrated-quantities-2-2" id="function-sum-integrated-quantities-2-2">function sum\_integrated\_quantities [2/2]</a>


```cpp
virtual void Nyx::sum_integrated_quantities () 
```



### <a href="#function-write-info-1-2" id="function-write-info-1-2">function write\_info [1/2]</a>


```cpp
void Nyx::write_info () 
```



### <a href="#function-write-info-1-2" id="function-write-info-1-2">function write\_info [1/2]</a>


```cpp
void Nyx::write_info () 
```


## Protected Static Functions Documentation


### <a href="#function-getcputime-1-2" id="function-getcputime-1-2">function getCPUTime [1/2]</a>


```cpp
static amrex::Real Nyx::getCPUTime () 
```



### <a href="#function-getcputime-2-2" id="function-getcputime-2-2">function getCPUTime [2/2]</a>


```cpp
static amrex::Real Nyx::getCPUTime () 
```



### <a href="#function-network-init-1-2" id="function-network-init-1-2">function network\_init [1/2]</a>


```cpp
static void Nyx::network_init () 
```



### <a href="#function-network-init-2-2" id="function-network-init-2-2">function network\_init [2/2]</a>


```cpp
static void Nyx::network_init () 
```



### <a href="#function-read-params-1-2" id="function-read-params-1-2">function read\_params [1/2]</a>


```cpp
static void Nyx::read_params () 
```



### <a href="#function-read-params-2-2" id="function-read-params-2-2">function read\_params [2/2]</a>


```cpp
static void Nyx::read_params () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Exec/Test_Only_Axions/Nyx.H`