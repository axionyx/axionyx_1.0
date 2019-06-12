
# Namespace network


[**Class List**](annotated.md) **>** [**network**](namespacenetwork.md)


















## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(kind=rt), dimension(nspec), save | [**aion**](namespacenetwork.md#variable-aion)  <br> |
|  real(kind=rt), dimension(nspec), save | [**ebin**](namespacenetwork.md#variable-ebin)  <br> |
|  integer, parameter | [**naux**](namespacenetwork.md#variable-naux)   = = 0<br> |
|  logical, save | [**network\_initialized**](namespacenetwork.md#variable-network-initialized)   = = .false.<br> |
|  integer, parameter | [**nspec**](namespacenetwork.md#variable-nspec)   = = 2<br> |
|  character(len=5), dimension(naux), save | [**short\_aux\_names**](namespacenetwork.md#variable-short-aux-names)  <br> |
|  character(len=5), dimension(nspec), save | [**short\_spec\_names**](namespacenetwork.md#variable-short-spec-names)  <br> |
|  character(len=16), dimension(nspec), save | [**spec\_names**](namespacenetwork.md#variable-spec-names)  <br> |
|  real(kind=rt), dimension(nspec), save | [**zion**](namespacenetwork.md#variable-zion)  <br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**network\_init**](namespacenetwork.md#function-network-init) () <br> |
|  integer function | [**network\_species\_index**](namespacenetwork.md#function-network-species-index) (name name) <br> |








## Public Attributes Documentation


### <a href="#variable-aion" id="variable-aion">variable aion </a>


```cpp
real(kind=rt), dimension(nspec), save network::aion;
```



### <a href="#variable-ebin" id="variable-ebin">variable ebin </a>


```cpp
real(kind=rt), dimension(nspec), save network::ebin;
```



### <a href="#variable-naux" id="variable-naux">variable naux </a>


```cpp
integer, parameter network::naux;
```



### <a href="#variable-network-initialized" id="variable-network-initialized">variable network\_initialized </a>


```cpp
logical, save network::network_initialized;
```



### <a href="#variable-nspec" id="variable-nspec">variable nspec </a>


```cpp
integer, parameter network::nspec;
```



### <a href="#variable-short-aux-names" id="variable-short-aux-names">variable short\_aux\_names </a>


```cpp
character (len= 5), dimension(naux), save network::short_aux_names;
```



### <a href="#variable-short-spec-names" id="variable-short-spec-names">variable short\_spec\_names </a>


```cpp
character (len= 5), dimension(nspec), save network::short_spec_names;
```



### <a href="#variable-spec-names" id="variable-spec-names">variable spec\_names </a>


```cpp
character (len=16), dimension(nspec), save network::spec_names;
```



### <a href="#variable-zion" id="variable-zion">variable zion </a>


```cpp
real(kind=rt), dimension(nspec), save network::zion;
```


## Public Functions Documentation


### <a href="#function-network-init" id="function-network-init">function network\_init </a>


```cpp
subroutine network::network_init () 
```



### <a href="#function-network-species-index" id="function-network-species-index">function network\_species\_index </a>


```cpp
integer function network::network_species_index (
    name name
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/Network/network.f90`