
# File distribution.c


[**File List**](files.md) **>** [**Exec**](dir_43a12cefb7942b6f49b5b628aafd3192.md) **>** [**Test\_Only\_Axions**](dir_eb24725df855cf6c732a19e4912f662a.md) **>** [**distribution.c**](distribution_8c.md)

[Go to the source code of this file.](distribution_8c_source.md)



* `#include <assert.h>`
* `#include <mpi.h>`
* `#include <stdbool.h>`
* `#include <stddef.h>`
* `#include <stdio.h>`
* `#include <stdlib.h>`
* `#include <stdint.h>`
* `#include "distribution_c.h"`











## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**@0**](distribution_8c.md#enum-@0)  <br> |




## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**Coord\_cube**](distribution_8c.md#function-coord-cube) (int myrank, int coord, distribution\_t \* d) <br> |
|  void | [**Coord\_x\_pencils**](distribution_8c.md#function-coord-x-pencils) (int myrank, int coord, distribution\_t \* d) <br> |
|  void | [**Coord\_y\_pencils**](distribution_8c.md#function-coord-y-pencils) (int myrank, int coord, distribution\_t \* d) <br> |
|  void | [**Coord\_z\_pencils**](distribution_8c.md#function-coord-z-pencils) (int myrank, int coord, distribution\_t \* d) <br> |
|  void | [**Custom3D\_Dims\_create**](distribution_8c.md#function-custom3d-dims-create) (const int Ndims, int nproc, int ndims, int dims) <br> |
|  void | [**Rank\_cube**](distribution_8c.md#function-rank-cube) (int \* myrank, int coord, distribution\_t \* d) <br> |
|  void | [**Rank\_x\_pencils**](distribution_8c.md#function-rank-x-pencils) (int \* myrank, int coord, distribution\_t \* d) <br> |
|  void | [**Rank\_y\_pencils**](distribution_8c.md#function-rank-y-pencils) (int \* myrank, int coord, distribution\_t \* d) <br> |
|  void | [**Rank\_z\_pencils**](distribution_8c.md#function-rank-z-pencils) (int \* myrank, int coord, distribution\_t \* d) <br> |
|  void | [**distribution\_1\_to\_3**](distribution_8c.md#function-distribution-1-to-3) (const complex\_t \* a, complex\_t \* b, distribution\_t \* d) <br> |
|  void | [**distribution\_2\_to\_3**](distribution_8c.md#function-distribution-2-to-3) (const complex\_t \* a, complex\_t \* b, distribution\_t \* d, int z\_dim) <br> |
|  void | [**distribution\_3\_to\_1**](distribution_8c.md#function-distribution-3-to-1) (const complex\_t \* a, complex\_t \* b, distribution\_t \* d) <br> |
|  void | [**distribution\_3\_to\_2**](distribution_8c.md#function-distribution-3-to-2) (const complex\_t \* a, complex\_t \* b, distribution\_t \* d, int z\_dim) <br> |
|  void | [**distribution\_assert\_commensurate**](distribution_8c.md#function-distribution-assert-commensurate) (distribution\_t \* d) <br> |
|  void | [**distribution\_fini**](distribution_8c.md#function-distribution-fini) (distribution\_t \* d) <br> |
|  void | [**distribution\_init**](distribution_8c.md#function-distribution-init) (MPI\_Comm comm, const int n, const int Ndims, distribution\_t \* d, const int \* rmap, bool debug) <br> |
|  void | [**distribution\_init\_explicit**](distribution_8c.md#function-distribution-init-explicit) (MPI\_Comm comm, const int n, int nproc\_1d, int nproc\_2d\_x, int nproc\_2d\_y, int nproc\_2d\_z, int nproc\_3d, distribution\_t \* d, bool debug) <br> |

## Public Static Functions

| Type | Name |
| ---: | :--- |
|  void | [**redistribute**](distribution_8c.md#function-redistribute) (const complex\_t \* a, complex\_t \* b, distribution\_t \* d, int direction) <br> |
|  void | [**redistribute\_2\_and\_3**](distribution_8c.md#function-redistribute-2-and-3) (const complex\_t \* a, complex\_t \* b, distribution\_t \* d, int direction, int z\_dim) <br> |
|  void | [**redistribute\_slab**](distribution_8c.md#function-redistribute-slab) (const complex\_t \* a, complex\_t \* b, distribution\_t \* d, int direction) <br> |
|  char \* | [**separator**](distribution_8c.md#function-separator) (int i, int n) <br> |






## Macros

| Type | Name |
| ---: | :--- |
| define  | [**BE\_VERBOSE**](distribution_8c.md#define-be-verbose)  () 0<br> |
| define  | [**DEBUG\_CONDITION**](distribution_8c.md#define-debug-condition)  () false<br> |
| define  | [**USE\_SLAB\_WORKAROUND**](distribution_8c.md#define-use-slab-workaround)  () 0<br> |

## Public Types Documentation


### <a href="#enum-@0" id="enum-@0">enum @0 </a>


```cpp
enum @0 {
    REDISTRIBUTE_1_TO_3,
    REDISTRIBUTE_3_TO_1,
    REDISTRIBUTE_2_TO_3,
    REDISTRIBUTE_3_TO_2
};
```


## Public Functions Documentation


### <a href="#function-coord-cube" id="function-coord-cube">function Coord\_cube </a>


```cpp
void Coord_cube (
    int myrank,
    int coord,
    distribution_t * d
) 
```



### <a href="#function-coord-x-pencils" id="function-coord-x-pencils">function Coord\_x\_pencils </a>


```cpp
void Coord_x_pencils (
    int myrank,
    int coord,
    distribution_t * d
) 
```



### <a href="#function-coord-y-pencils" id="function-coord-y-pencils">function Coord\_y\_pencils </a>


```cpp
void Coord_y_pencils (
    int myrank,
    int coord,
    distribution_t * d
) 
```



### <a href="#function-coord-z-pencils" id="function-coord-z-pencils">function Coord\_z\_pencils </a>


```cpp
void Coord_z_pencils (
    int myrank,
    int coord,
    distribution_t * d
) 
```



### <a href="#function-custom3d-dims-create" id="function-custom3d-dims-create">function Custom3D\_Dims\_create </a>


```cpp
void Custom3D_Dims_create (
    const int Ndims,
    int nproc,
    int ndims,
    int dims
) 
```



### <a href="#function-rank-cube" id="function-rank-cube">function Rank\_cube </a>


```cpp
void Rank_cube (
    int * myrank,
    int coord,
    distribution_t * d
) 
```



### <a href="#function-rank-x-pencils" id="function-rank-x-pencils">function Rank\_x\_pencils </a>


```cpp
void Rank_x_pencils (
    int * myrank,
    int coord,
    distribution_t * d
) 
```



### <a href="#function-rank-y-pencils" id="function-rank-y-pencils">function Rank\_y\_pencils </a>


```cpp
void Rank_y_pencils (
    int * myrank,
    int coord,
    distribution_t * d
) 
```



### <a href="#function-rank-z-pencils" id="function-rank-z-pencils">function Rank\_z\_pencils </a>


```cpp
void Rank_z_pencils (
    int * myrank,
    int coord,
    distribution_t * d
) 
```



### <a href="#function-distribution-1-to-3" id="function-distribution-1-to-3">function distribution\_1\_to\_3 </a>


```cpp
void distribution_1_to_3 (
    const complex_t * a,
    complex_t * b,
    distribution_t * d
) 
```



### <a href="#function-distribution-2-to-3" id="function-distribution-2-to-3">function distribution\_2\_to\_3 </a>


```cpp
void distribution_2_to_3 (
    const complex_t * a,
    complex_t * b,
    distribution_t * d,
    int z_dim
) 
```



### <a href="#function-distribution-3-to-1" id="function-distribution-3-to-1">function distribution\_3\_to\_1 </a>


```cpp
void distribution_3_to_1 (
    const complex_t * a,
    complex_t * b,
    distribution_t * d
) 
```



### <a href="#function-distribution-3-to-2" id="function-distribution-3-to-2">function distribution\_3\_to\_2 </a>


```cpp
void distribution_3_to_2 (
    const complex_t * a,
    complex_t * b,
    distribution_t * d,
    int z_dim
) 
```



### <a href="#function-distribution-assert-commensurate" id="function-distribution-assert-commensurate">function distribution\_assert\_commensurate </a>


```cpp
void distribution_assert_commensurate (
    distribution_t * d
) 
```



### <a href="#function-distribution-fini" id="function-distribution-fini">function distribution\_fini </a>


```cpp
void distribution_fini (
    distribution_t * d
) 
```



### <a href="#function-distribution-init" id="function-distribution-init">function distribution\_init </a>


```cpp
void distribution_init (
    MPI_Comm comm,
    const int n,
    const int Ndims,
    distribution_t * d,
    const int * rmap,
    bool debug
) 
```



### <a href="#function-distribution-init-explicit" id="function-distribution-init-explicit">function distribution\_init\_explicit </a>


```cpp
void distribution_init_explicit (
    MPI_Comm comm,
    const int n,
    int nproc_1d,
    int nproc_2d_x,
    int nproc_2d_y,
    int nproc_2d_z,
    int nproc_3d,
    distribution_t * d,
    bool debug
) 
```


## Public Static Functions Documentation


### <a href="#function-redistribute" id="function-redistribute">function redistribute </a>


```cpp
static void redistribute (
    const complex_t * a,
    complex_t * b,
    distribution_t * d,
    int direction
) 
```



### <a href="#function-redistribute-2-and-3" id="function-redistribute-2-and-3">function redistribute\_2\_and\_3 </a>


```cpp
static void redistribute_2_and_3 (
    const complex_t * a,
    complex_t * b,
    distribution_t * d,
    int direction,
    int z_dim
) 
```



### <a href="#function-redistribute-slab" id="function-redistribute-slab">function redistribute\_slab </a>


```cpp
static void redistribute_slab (
    const complex_t * a,
    complex_t * b,
    distribution_t * d,
    int direction
) 
```



### <a href="#function-separator" id="function-separator">function separator </a>


```cpp
static inline char * separator (
    int i,
    int n
) 
```

## Macro Definition Documentation



### <a href="#define-be-verbose" id="define-be-verbose">define BE\_VERBOSE </a>


```cpp
#define BE_VERBOSE () 
```



### <a href="#define-debug-condition" id="define-debug-condition">define DEBUG\_CONDITION </a>


```cpp
#define DEBUG_CONDITION () 
```



### <a href="#define-use-slab-workaround" id="define-use-slab-workaround">define USE\_SLAB\_WORKAROUND </a>


```cpp
#define USE_SLAB_WORKAROUND () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Exec/Test_Only_Axions/distribution.c`