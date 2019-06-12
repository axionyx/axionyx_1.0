
# Class amrex::FillPatchIteratorHelper


[**Class List**](annotated.md) **>** [**amrex**](namespaceamrex.md) **>** [**FillPatchIteratorHelper**](classamrex_1_1FillPatchIteratorHelper.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**FillPatchIteratorHelper**](classamrex_1_1FillPatchIteratorHelper.md#function-fillpatchiteratorhelper-1-4) ([**AmrLevel**](classamrex_1_1AmrLevel.md) & amrlevel, MultiFab & leveldata) <br> |
|   | [**FillPatchIteratorHelper**](classamrex_1_1FillPatchIteratorHelper.md#function-fillpatchiteratorhelper-2-4) ([**AmrLevel**](classamrex_1_1AmrLevel.md) & amrlevel, MultiFab & leveldata, int boxGrow, Real time, int state\_indx, int scomp, int ncomp, Interpolater \* mapper) <br> |
|  void | [**Initialize**](classamrex_1_1FillPatchIteratorHelper.md#function-initialize) (int boxGrow, Real time, int state\_indx, int scomp, int ncomp, Interpolater \* mapper) <br> |
|  void | [**fill**](classamrex_1_1FillPatchIteratorHelper.md#function-fill) (FArrayBox & fab, int dcomp, int idx) <br> |
|   | [**~FillPatchIteratorHelper**](classamrex_1_1FillPatchIteratorHelper.md#function-fillpatchiteratorhelper) () <br> |








## Public Functions Documentation


### <a href="#function-fillpatchiteratorhelper-1-4" id="function-fillpatchiteratorhelper-1-4">function FillPatchIteratorHelper [1/4]</a>


```cpp
amrex::FillPatchIteratorHelper::FillPatchIteratorHelper (
    AmrLevel & amrlevel,
    MultiFab & leveldata
) 
```



### <a href="#function-fillpatchiteratorhelper-2-4" id="function-fillpatchiteratorhelper-2-4">function FillPatchIteratorHelper [2/4]</a>


```cpp
amrex::FillPatchIteratorHelper::FillPatchIteratorHelper (
    AmrLevel & amrlevel,
    MultiFab & leveldata,
    int boxGrow,
    Real time,
    int state_indx,
    int scomp,
    int ncomp,
    Interpolater * mapper
) 
```



### <a href="#function-initialize" id="function-initialize">function Initialize </a>


```cpp
void amrex::FillPatchIteratorHelper::Initialize (
    int boxGrow,
    Real time,
    int state_indx,
    int scomp,
    int ncomp,
    Interpolater * mapper
) 
```



### <a href="#function-fill" id="function-fill">function fill </a>


```cpp
void amrex::FillPatchIteratorHelper::fill (
    FArrayBox & fab,
    int dcomp,
    int idx
) 
```



### <a href="#function-fillpatchiteratorhelper" id="function-fillpatchiteratorhelper">function ~FillPatchIteratorHelper </a>


```cpp
amrex::FillPatchIteratorHelper::~FillPatchIteratorHelper () 
```

## Friends Documentation



### <a href="#friend-fillpatchiterator" id="friend-fillpatchiterator">friend FillPatchIterator </a>


```cpp
friend class amrex::FillPatchIteratorHelper::FillPatchIterator () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/AMReX_axionyx/AMReX_AmrLevel.H`