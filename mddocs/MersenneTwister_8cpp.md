
# File MersenneTwister.cpp


[**File List**](files.md) **>** [**Forcing**](dir_45682215f16eaf57f766b3c547de68bc.md) **>** [**MersenneTwister.cpp**](MersenneTwister_8cpp.md)

[Go to the source code of this file.](MersenneTwister_8cpp_source.md)



* `#include <stdio.h>`
* `#include <stdlib.h>`
* `#include <math.h>`
* `#include <iostream>`
* `#include <fstream>`













## Public Attributes

| Type | Name |
| ---: | :--- |
|  unsigned long int | [**mt\_buffer**](MersenneTwister_8cpp.md#variable-mt-buffer)  <br> |
|  int | [**mt\_index**](MersenneTwister_8cpp.md#variable-mt-index)  <br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**mt\_init**](MersenneTwister_8cpp.md#function-mt-init) (unsigned int seed) <br> |
|  unsigned long int | [**mt\_random**](MersenneTwister_8cpp.md#function-mt-random) () <br> |
|  void | [**mt\_read**](MersenneTwister_8cpp.md#function-mt-read) (std::ifstream & input) <br> |
|  void | [**mt\_write**](MersenneTwister_8cpp.md#function-mt-write) (std::ofstream & output) <br> |







## Macros

| Type | Name |
| ---: | :--- |
| define  | [**LOWER\_MASK**](MersenneTwister_8cpp.md#define-lower-mask)  () 0x7FFFFFFF<br> |
| define  | [**MAGIC**](MersenneTwister_8cpp.md#define-magic) (s) (s) (((s)&1)\*MATRIX\_A)<br> |
| define  | [**MATRIX\_A**](MersenneTwister_8cpp.md#define-matrix-a)  () 0x9908B0DF<br> |
| define  | [**MT\_IA**](MersenneTwister_8cpp.md#define-mt-ia)  () 397<br> |
| define  | [**MT\_IB**](MersenneTwister_8cpp.md#define-mt-ib)  () (MT\_LEN - MT\_IA)<br> |
| define  | [**MT\_LEN**](MersenneTwister_8cpp.md#define-mt-len)  () 624<br> |
| define  | [**TWIST**](MersenneTwister_8cpp.md#define-twist) (b, i, j) (b, i, j) ((b)[i] & UPPER\_MASK) | ((b)[j] & LOWER\_MASK)<br> |
| define  | [**UPPER\_MASK**](MersenneTwister_8cpp.md#define-upper-mask)  () 0x80000000<br> |

## Public Attributes Documentation


### <a href="#variable-mt-buffer" id="variable-mt-buffer">variable mt\_buffer </a>


```cpp
unsigned long int mt_buffer[MT_LEN];
```



### <a href="#variable-mt-index" id="variable-mt-index">variable mt\_index </a>


```cpp
int mt_index;
```


## Public Functions Documentation


### <a href="#function-mt-init" id="function-mt-init">function mt\_init </a>


```cpp
void mt_init (
    unsigned int seed
) 
```



### <a href="#function-mt-random" id="function-mt-random">function mt\_random </a>


```cpp
unsigned long int mt_random () 
```



### <a href="#function-mt-read" id="function-mt-read">function mt\_read </a>


```cpp
void mt_read (
    std::ifstream & input
) 
```



### <a href="#function-mt-write" id="function-mt-write">function mt\_write </a>


```cpp
void mt_write (
    std::ofstream & output
) 
```

## Macro Definition Documentation



### <a href="#define-lower-mask" id="define-lower-mask">define LOWER\_MASK </a>


```cpp
#define LOWER_MASK () 
```



### <a href="#define-magic" id="define-magic">define MAGIC </a>


```cpp
#define MAGIC (
    s
) (((s)&1)*MATRIX_A)
```



### <a href="#define-matrix-a" id="define-matrix-a">define MATRIX\_A </a>


```cpp
#define MATRIX_A () 
```



### <a href="#define-mt-ia" id="define-mt-ia">define MT\_IA </a>


```cpp
#define MT_IA () 
```



### <a href="#define-mt-ib" id="define-mt-ib">define MT\_IB </a>


```cpp
#define MT_IB () 
```



### <a href="#define-mt-len" id="define-mt-len">define MT\_LEN </a>


```cpp
#define MT_LEN () 
```



### <a href="#define-twist" id="define-twist">define TWIST </a>


```cpp
#define TWIST (
    b,
    i,
    j
) ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
```



### <a href="#define-upper-mask" id="define-upper-mask">define UPPER\_MASK </a>


```cpp
#define UPPER_MASK () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/Forcing/MersenneTwister.cpp`