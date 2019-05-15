
# Namespace forcing\_spect\_module


[**Class List**](annotated.md) **>** [**forcing\_spect\_module**](namespaceforcing__spect__module.md)


















## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter | [**lenvec**](namespaceforcing__spect__module.md#variable-lenvec)   = = 4<br> |
|  real(rt), dimension(:,:), allocatable | [**modes\_even**](namespaceforcing__spect__module.md#variable-modes-even)  <br> |
|  real(rt), dimension(:,:), allocatable | [**modes\_odd**](namespaceforcing__spect__module.md#variable-modes-odd)  <br> |
|  integer | [**num\_modes**](namespaceforcing__spect__module.md#variable-num-modes)  <br> |
|  integer | [**num\_modes\_ext**](namespaceforcing__spect__module.md#variable-num-modes-ext)  <br> |
|  integer, dimension(:,:), allocatable | [**wavevectors**](namespaceforcing__spect__module.md#variable-wavevectors)  <br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**fort\_alloc\_spect**](namespaceforcing__spect__module.md#function-fort-alloc-spect) (length length) <br> |
|  subroutine, public | [**fort\_set\_modes**](namespaceforcing__spect__module.md#function-fort-set-modes) (even even, odd odd, length length, comp comp) <br> |
|  subroutine, public | [**fort\_set\_wavevector**](namespaceforcing__spect__module.md#function-fort-set-wavevector) (kvect kvect, m m) <br> |








## Public Attributes Documentation


### <a href="#variable-lenvec" id="variable-lenvec">variable lenvec </a>


```cpp
integer, parameter forcing_spect_module::lenvec;
```



### <a href="#variable-modes-even" id="variable-modes-even">variable modes\_even </a>


```cpp
real(rt), dimension(:,:), allocatable forcing_spect_module::modes_even;
```



### <a href="#variable-modes-odd" id="variable-modes-odd">variable modes\_odd </a>


```cpp
real(rt), dimension(:,:), allocatable forcing_spect_module::modes_odd;
```



### <a href="#variable-num-modes" id="variable-num-modes">variable num\_modes </a>


```cpp
integer forcing_spect_module::num_modes;
```



### <a href="#variable-num-modes-ext" id="variable-num-modes-ext">variable num\_modes\_ext </a>


```cpp
integer forcing_spect_module::num_modes_ext;
```



### <a href="#variable-wavevectors" id="variable-wavevectors">variable wavevectors </a>


```cpp
integer, dimension(:,:), allocatable forcing_spect_module::wavevectors;
```


## Public Functions Documentation


### <a href="#function-fort-alloc-spect" id="function-fort-alloc-spect">function fort\_alloc\_spect </a>


```cpp
subroutine, public forcing_spect_module::fort_alloc_spect (
    length length
) 
```



### <a href="#function-fort-set-modes" id="function-fort-set-modes">function fort\_set\_modes </a>


```cpp
subroutine, public forcing_spect_module::fort_set_modes (
    even even,
    odd odd,
    length length,
    comp comp
) 
```



### <a href="#function-fort-set-wavevector" id="function-fort-set-wavevector">function fort\_set\_wavevector </a>


```cpp
subroutine, public forcing_spect_module::fort_set_wavevector (
    kvect kvect,
    m m
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/Forcing/forcing_spect.f90`