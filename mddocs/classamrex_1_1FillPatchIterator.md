
# Class amrex::FillPatchIterator


[**Class List**](annotated.md) **>** [**amrex**](namespaceamrex.md) **>** [**FillPatchIterator**](classamrex_1_1FillPatchIterator.md)








Inherits the following classes: MFIter












## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**FillPatchIterator**](classamrex_1_1FillPatchIterator.md#function-fillpatchiterator-1-4) ([**AmrLevel**](classamrex_1_1AmrLevel.md) & amrlevel, MultiFab & leveldata) <br> |
|   | [**FillPatchIterator**](classamrex_1_1FillPatchIterator.md#function-fillpatchiterator-2-4) ([**AmrLevel**](classamrex_1_1AmrLevel.md) & amrlevel, MultiFab & leveldata, int boxGrow, Real time, int state\_indx, int scomp, int ncomp) <br> |
|  void | [**Initialize**](classamrex_1_1FillPatchIterator.md#function-initialize) (int boxGrow, Real time, int state\_indx, int scomp, int ncomp) <br> |
|  Box | [**UngrownBox**](classamrex_1_1FillPatchIterator.md#function-ungrownbox) () noexcept const<br> |
|  MultiFab & | [**get\_mf**](classamrex_1_1FillPatchIterator.md#function-get-mf) () noexcept<br> |
|  FArrayBox & | [**operator()**](classamrex_1_1FillPatchIterator.md#function-operator()) () noexcept<br> |
|   | [**~FillPatchIterator**](classamrex_1_1FillPatchIterator.md#function-fillpatchiterator) () <br> |








## Public Functions Documentation


### <a href="#function-fillpatchiterator-1-4" id="function-fillpatchiterator-1-4">function FillPatchIterator [1/4]</a>


```cpp
amrex::FillPatchIterator::FillPatchIterator (
    AmrLevel & amrlevel,
    MultiFab & leveldata
) 
```



### <a href="#function-fillpatchiterator-2-4" id="function-fillpatchiterator-2-4">function FillPatchIterator [2/4]</a>


```cpp
amrex::FillPatchIterator::FillPatchIterator (
    AmrLevel & amrlevel,
    MultiFab & leveldata,
    int boxGrow,
    Real time,
    int state_indx,
    int scomp,
    int ncomp
) 
```



### <a href="#function-initialize" id="function-initialize">function Initialize </a>


```cpp
void amrex::FillPatchIterator::Initialize (
    int boxGrow,
    Real time,
    int state_indx,
    int scomp,
    int ncomp
) 
```



### <a href="#function-ungrownbox" id="function-ungrownbox">function UngrownBox </a>


```cpp
inline Box amrex::FillPatchIterator::UngrownBox () noexcept const
```



### <a href="#function-get-mf" id="function-get-mf">function get\_mf </a>


```cpp
inline MultiFab & amrex::FillPatchIterator::get_mf () noexcept
```



### <a href="#function-operator()" id="function-operator()">function operator() </a>


```cpp
inline FArrayBox & amrex::FillPatchIterator::operator() () noexcept
```



### <a href="#function-fillpatchiterator" id="function-fillpatchiterator">function ~FillPatchIterator </a>


```cpp
amrex::FillPatchIterator::~FillPatchIterator () 
```

## Friends Documentation



### <a href="#friend-amrlevel" id="friend-amrlevel">friend AmrLevel </a>


```cpp
friend class amrex::FillPatchIterator::AmrLevel () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/AMReX_axionyx/AMReX_AmrLevel.H`