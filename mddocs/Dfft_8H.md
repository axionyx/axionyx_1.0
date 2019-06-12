
# File Dfft.H


[**File List**](files.md) **>** [**Exec**](dir_43a12cefb7942b6f49b5b628aafd3192.md) **>** [**Test\_Only\_Axions**](dir_eb24725df855cf6c732a19e4912f662a.md) **>** [**Dfft.H**](Dfft_8H.md)

[Go to the source code of this file.](Dfft_8H_source.md)



* `#include <fftw3.h>`
* `#include "complex-type.h"`
* `#include "Distribution.H"`
* `#include "Error.h"`









## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**hacc**](namespacehacc.md) <br> |

## Classes

| Type | Name |
| ---: | :--- |
| class | [**Dfft**](classhacc_1_1Dfft.md) <br> |












## Macros

| Type | Name |
| ---: | :--- |
| define  | [**DFFT\_TIMING**](Dfft_8H.md#define-dfft-timing)  () 0<br> |
| define  | [**FFTW\_ADDR**](Dfft_8H.md#define-fftw-addr) (X) (X) reinterpret\_cast&lt;fftw\_complex\*&gt;(&(X)[0])<br> |
## Macro Definition Documentation



### <a href="#define-dfft-timing" id="define-dfft-timing">define DFFT\_TIMING </a>


```cpp
#define DFFT_TIMING () 
```



### <a href="#define-fftw-addr" id="define-fftw-addr">define FFTW\_ADDR </a>


```cpp
#define FFTW_ADDR (
    X
) reinterpret_cast<fftw_complex*>(&(X)[0])
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Exec/Test_Only_Axions/Dfft.H`