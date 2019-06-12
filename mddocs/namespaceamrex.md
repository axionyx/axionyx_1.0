
# Namespace amrex


[**Class List**](annotated.md) **>** [**amrex**](namespaceamrex.md)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**anonymous\_namespace{AMReX\_Amr.cpp}**](namespaceamrex_1_1anonymous__namespace_02AMReX__Amr_8cpp_03.md) <br> |
| namespace | [**anonymous\_namespace{AMReX\_AmrCore.cpp}**](namespaceamrex_1_1anonymous__namespace_02AMReX__AmrCore_8cpp_03.md) <br> |

## Classes

| Type | Name |
| ---: | :--- |
| class | [**Amr**](classamrex_1_1Amr.md) <br>_Manage hierarchy of levels for time-dependent AMR computations._  |
| class | [**AmrCore**](classamrex_1_1AmrCore.md) <br>_Provide basic functionalities to set up an AMR hierarchy._  |
| class | [**AmrLevel**](classamrex_1_1AmrLevel.md) <br>_Virtual base class for managing individual levels._ [_**AmrLevel**_](classamrex_1_1AmrLevel.md) _functions both as a container for state data on a level and also manages the advancement of data in time._ |
| class | [**FillPatchIterator**](classamrex_1_1FillPatchIterator.md) <br> |
| class | [**FillPatchIteratorHelper**](classamrex_1_1FillPatchIteratorHelper.md) <br> |
| class | [**MFGraph**](classamrex_1_1MFGraph.md) &lt;class T&gt;<br> |



## Public Attributes

| Type | Name |
| ---: | :--- |
|  bool | [**checkpoint\_files\_output**](namespaceamrex.md#variable-checkpoint-files-output)  <br> |
|  int | [**checkpoint\_nfiles**](namespaceamrex.md#variable-checkpoint-nfiles)  <br> |
|  int | [**checkpoint\_on\_restart**](namespaceamrex.md#variable-checkpoint-on-restart)  <br> |
|  int | [**compute\_new\_dt\_on\_regrid**](namespaceamrex.md#variable-compute-new-dt-on-regrid)  <br> |
|  int | [**insitu\_on\_restart**](namespaceamrex.md#variable-insitu-on-restart)  <br> |
|  int | [**mffile\_nstreams**](namespaceamrex.md#variable-mffile-nstreams)  <br> |
|  bool | [**plot\_files\_output**](namespaceamrex.md#variable-plot-files-output)  <br> |
|  int | [**plot\_nfiles**](namespaceamrex.md#variable-plot-nfiles)  <br> |
|  int | [**plotfile\_on\_restart**](namespaceamrex.md#variable-plotfile-on-restart)  <br> |
|  bool | [**precreateDirectories**](namespaceamrex.md#variable-precreatedirectories)  <br> |
|  bool | [**prereadFAHeaders**](namespaceamrex.md#variable-prereadfaheaders)  <br> |
|  int | [**probinit\_natonce**](namespaceamrex.md#variable-probinit-natonce)  <br> |
|  int | [**regrid\_on\_restart**](namespaceamrex.md#variable-regrid-on-restart)  <br> |
|  int | [**use\_efficient\_regrid**](namespaceamrex.md#variable-use-efficient-regrid)  <br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|  const char \* | [**buildInfoGetAMReXDir**](namespaceamrex.md#function-buildinfogetamrexdir) () <br> |
|  const char \* | [**buildInfoGetAux**](namespaceamrex.md#function-buildinfogetaux) (int i) <br> |
|  const char \* | [**buildInfoGetBuildDate**](namespaceamrex.md#function-buildinfogetbuilddate) () <br> |
|  const char \* | [**buildInfoGetBuildDir**](namespaceamrex.md#function-buildinfogetbuilddir) () <br> |
|  const char \* | [**buildInfoGetBuildGitHash**](namespaceamrex.md#function-buildinfogetbuildgithash) () <br> |
|  const char \* | [**buildInfoGetBuildGitName**](namespaceamrex.md#function-buildinfogetbuildgitname) () <br> |
|  const char \* | [**buildInfoGetBuildMachine**](namespaceamrex.md#function-buildinfogetbuildmachine) () <br> |
|  const char \* | [**buildInfoGetCXXFlags**](namespaceamrex.md#function-buildinfogetcxxflags) () <br> |
|  const char \* | [**buildInfoGetCXXName**](namespaceamrex.md#function-buildinfogetcxxname) () <br> |
|  const char \* | [**buildInfoGetComp**](namespaceamrex.md#function-buildinfogetcomp) () <br> |
|  const char \* | [**buildInfoGetCompVersion**](namespaceamrex.md#function-buildinfogetcompversion) () <br> |
|  const char \* | [**buildInfoGetFFlags**](namespaceamrex.md#function-buildinfogetfflags) () <br> |
|  const char \* | [**buildInfoGetFName**](namespaceamrex.md#function-buildinfogetfname) () <br> |
|  const char \* | [**buildInfoGetFcomp**](namespaceamrex.md#function-buildinfogetfcomp) () <br> |
|  const char \* | [**buildInfoGetFcompVersion**](namespaceamrex.md#function-buildinfogetfcompversion) () <br> |
|  const char \* | [**buildInfoGetGitHash**](namespaceamrex.md#function-buildinfogetgithash) (int i) <br> |
|  const char \* | [**buildInfoGetLibraries**](namespaceamrex.md#function-buildinfogetlibraries) () <br> |
|  const char \* | [**buildInfoGetLinkFlags**](namespaceamrex.md#function-buildinfogetlinkflags) () <br> |
|  const char \* | [**buildInfoGetModuleName**](namespaceamrex.md#function-buildinfogetmodulename) (int i) <br> |
|  const char \* | [**buildInfoGetModuleVal**](namespaceamrex.md#function-buildinfogetmoduleval) (int i) <br> |
|  int | [**buildInfoGetNumModules**](namespaceamrex.md#function-buildinfogetnummodules) () <br> |
|  VisMF::Header::Version | [**checkpoint\_headerversion**](namespaceamrex.md#function-checkpoint-headerversion) (VisMF::Header::Version\_v1) <br> |
|  VisMF::Header::Version | [**plot\_headerversion**](namespaceamrex.md#function-plot-headerversion) (VisMF::Header::Version\_v1) <br> |








## Public Attributes Documentation


### <a href="#variable-checkpoint-files-output" id="variable-checkpoint-files-output">variable checkpoint\_files\_output </a>


```cpp
bool amrex::checkpoint_files_output;
```



### <a href="#variable-checkpoint-nfiles" id="variable-checkpoint-nfiles">variable checkpoint\_nfiles </a>


```cpp
int amrex::checkpoint_nfiles;
```



### <a href="#variable-checkpoint-on-restart" id="variable-checkpoint-on-restart">variable checkpoint\_on\_restart </a>


```cpp
int amrex::checkpoint_on_restart;
```



### <a href="#variable-compute-new-dt-on-regrid" id="variable-compute-new-dt-on-regrid">variable compute\_new\_dt\_on\_regrid </a>


```cpp
int amrex::compute_new_dt_on_regrid;
```



### <a href="#variable-insitu-on-restart" id="variable-insitu-on-restart">variable insitu\_on\_restart </a>


```cpp
int amrex::insitu_on_restart;
```



### <a href="#variable-mffile-nstreams" id="variable-mffile-nstreams">variable mffile\_nstreams </a>


```cpp
int amrex::mffile_nstreams;
```



### <a href="#variable-plot-files-output" id="variable-plot-files-output">variable plot\_files\_output </a>


```cpp
bool amrex::plot_files_output;
```



### <a href="#variable-plot-nfiles" id="variable-plot-nfiles">variable plot\_nfiles </a>


```cpp
int amrex::plot_nfiles;
```



### <a href="#variable-plotfile-on-restart" id="variable-plotfile-on-restart">variable plotfile\_on\_restart </a>


```cpp
int amrex::plotfile_on_restart;
```



### <a href="#variable-precreatedirectories" id="variable-precreatedirectories">variable precreateDirectories </a>


```cpp
bool amrex::precreateDirectories;
```



### <a href="#variable-prereadfaheaders" id="variable-prereadfaheaders">variable prereadFAHeaders </a>


```cpp
bool amrex::prereadFAHeaders;
```



### <a href="#variable-probinit-natonce" id="variable-probinit-natonce">variable probinit\_natonce </a>


```cpp
int amrex::probinit_natonce;
```



### <a href="#variable-regrid-on-restart" id="variable-regrid-on-restart">variable regrid\_on\_restart </a>


```cpp
int amrex::regrid_on_restart;
```



### <a href="#variable-use-efficient-regrid" id="variable-use-efficient-regrid">variable use\_efficient\_regrid </a>


```cpp
int amrex::use_efficient_regrid;
```


## Public Functions Documentation


### <a href="#function-buildinfogetamrexdir" id="function-buildinfogetamrexdir">function buildInfoGetAMReXDir </a>


```cpp
const char * amrex::buildInfoGetAMReXDir () 
```



### <a href="#function-buildinfogetaux" id="function-buildinfogetaux">function buildInfoGetAux </a>


```cpp
const char * amrex::buildInfoGetAux (
    int i
) 
```



### <a href="#function-buildinfogetbuilddate" id="function-buildinfogetbuilddate">function buildInfoGetBuildDate </a>


```cpp
const char * amrex::buildInfoGetBuildDate () 
```



### <a href="#function-buildinfogetbuilddir" id="function-buildinfogetbuilddir">function buildInfoGetBuildDir </a>


```cpp
const char * amrex::buildInfoGetBuildDir () 
```



### <a href="#function-buildinfogetbuildgithash" id="function-buildinfogetbuildgithash">function buildInfoGetBuildGitHash </a>


```cpp
const char * amrex::buildInfoGetBuildGitHash () 
```



### <a href="#function-buildinfogetbuildgitname" id="function-buildinfogetbuildgitname">function buildInfoGetBuildGitName </a>


```cpp
const char * amrex::buildInfoGetBuildGitName () 
```



### <a href="#function-buildinfogetbuildmachine" id="function-buildinfogetbuildmachine">function buildInfoGetBuildMachine </a>


```cpp
const char * amrex::buildInfoGetBuildMachine () 
```



### <a href="#function-buildinfogetcxxflags" id="function-buildinfogetcxxflags">function buildInfoGetCXXFlags </a>


```cpp
const char * amrex::buildInfoGetCXXFlags () 
```



### <a href="#function-buildinfogetcxxname" id="function-buildinfogetcxxname">function buildInfoGetCXXName </a>


```cpp
const char * amrex::buildInfoGetCXXName () 
```



### <a href="#function-buildinfogetcomp" id="function-buildinfogetcomp">function buildInfoGetComp </a>


```cpp
const char * amrex::buildInfoGetComp () 
```



### <a href="#function-buildinfogetcompversion" id="function-buildinfogetcompversion">function buildInfoGetCompVersion </a>


```cpp
const char * amrex::buildInfoGetCompVersion () 
```



### <a href="#function-buildinfogetfflags" id="function-buildinfogetfflags">function buildInfoGetFFlags </a>


```cpp
const char * amrex::buildInfoGetFFlags () 
```



### <a href="#function-buildinfogetfname" id="function-buildinfogetfname">function buildInfoGetFName </a>


```cpp
const char * amrex::buildInfoGetFName () 
```



### <a href="#function-buildinfogetfcomp" id="function-buildinfogetfcomp">function buildInfoGetFcomp </a>


```cpp
const char * amrex::buildInfoGetFcomp () 
```



### <a href="#function-buildinfogetfcompversion" id="function-buildinfogetfcompversion">function buildInfoGetFcompVersion </a>


```cpp
const char * amrex::buildInfoGetFcompVersion () 
```



### <a href="#function-buildinfogetgithash" id="function-buildinfogetgithash">function buildInfoGetGitHash </a>


```cpp
const char * amrex::buildInfoGetGitHash (
    int i
) 
```



### <a href="#function-buildinfogetlibraries" id="function-buildinfogetlibraries">function buildInfoGetLibraries </a>


```cpp
const char * amrex::buildInfoGetLibraries () 
```



### <a href="#function-buildinfogetlinkflags" id="function-buildinfogetlinkflags">function buildInfoGetLinkFlags </a>


```cpp
const char * amrex::buildInfoGetLinkFlags () 
```



### <a href="#function-buildinfogetmodulename" id="function-buildinfogetmodulename">function buildInfoGetModuleName </a>


```cpp
const char * amrex::buildInfoGetModuleName (
    int i
) 
```



### <a href="#function-buildinfogetmoduleval" id="function-buildinfogetmoduleval">function buildInfoGetModuleVal </a>


```cpp
const char * amrex::buildInfoGetModuleVal (
    int i
) 
```



### <a href="#function-buildinfogetnummodules" id="function-buildinfogetnummodules">function buildInfoGetNumModules </a>


```cpp
int amrex::buildInfoGetNumModules () 
```



### <a href="#function-checkpoint-headerversion" id="function-checkpoint-headerversion">function checkpoint\_headerversion </a>


```cpp
VisMF::Header::Version amrex::checkpoint_headerversion (
    VisMF::Header::Version_v1
) 
```



### <a href="#function-plot-headerversion" id="function-plot-headerversion">function plot\_headerversion </a>


```cpp
VisMF::Header::Version amrex::plot_headerversion (
    VisMF::Header::Version_v1
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Exec/Test_Only_Axions/AMReX_buildInfo.cpp`