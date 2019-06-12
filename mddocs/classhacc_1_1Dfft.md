
# Class hacc::Dfft


[**Class List**](annotated.md) **>** [**hacc**](namespacehacc.md) **>** [**Dfft**](classhacc_1_1Dfft.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Dfft**](classhacc_1_1Dfft.md#function-dfft-1-2) (Distribution & dist) <br> |
|   | [**Dfft**](classhacc_1_1Dfft.md#function-dfft-2-2) (Distribution & dist, complex\_t \* forward\_output, complex\_t \* forward\_scratch, complex\_t \* backward\_input, complex\_t \* backward\_scratch, unsigned int flags=FFTW\_MEASURE) <br> |
|  void | [**backward**](classhacc_1_1Dfft.md#function-backward-1-3) (complex\_t \* out) <br> |
|  void | [**backward**](classhacc_1_1Dfft.md#function-backward-2-3) (float \* out, size\_t ghost0, size\_t ghost1) <br> |
|  void | [**backward**](classhacc_1_1Dfft.md#function-backward-3-3) (float \* out, size\_t ghost) <br> |
|  MPI\_Comm | [**cartcomm\_kspace**](classhacc_1_1Dfft.md#function-cartcomm-kspace) () const<br> |
|  MPI\_Comm | [**cartcomm\_rspace**](classhacc_1_1Dfft.md#function-cartcomm-rspace) () const<br> |
|  void | [**forward**](classhacc_1_1Dfft.md#function-forward-1-3) (complex\_t const \* in) <br> |
|  void | [**forward**](classhacc_1_1Dfft.md#function-forward-2-3) (float const \* in, size\_t ghost0, size\_t ghost1) <br> |
|  void | [**forward**](classhacc_1_1Dfft.md#function-forward-3-3) (float const \* in, size\_t ghost) <br> |
|  Distribution & | [**get\_d**](classhacc_1_1Dfft.md#function-get-d) () <br> |
|  int | [**global\_ng**](classhacc_1_1Dfft.md#function-global-ng-1-2) (int i) const<br> |
|  int const (& | [**global\_ng**](classhacc_1_1Dfft.md#function-global-ng-2-2) () <br> |
|  size\_t | [**global\_size**](classhacc_1_1Dfft.md#function-global-size) () const<br> |
|  int | [**local\_ng\_kspace**](classhacc_1_1Dfft.md#function-local-ng-kspace-1-2) (int i) const<br> |
|  int const (& | [**local\_ng\_kspace**](classhacc_1_1Dfft.md#function-local-ng-kspace-2-2) () <br> |
|  int | [**local\_ng\_rspace**](classhacc_1_1Dfft.md#function-local-ng-rspace-1-2) (int i) const<br> |
|  int const (& | [**local\_ng\_rspace**](classhacc_1_1Dfft.md#function-local-ng-rspace-2-2) () <br> |
|  size\_t | [**local\_size**](classhacc_1_1Dfft.md#function-local-size) () const<br> |
|  void | [**makePlans**](classhacc_1_1Dfft.md#function-makeplans) (complex\_t \* forward\_output, complex\_t \* forward\_scratch, complex\_t \* backward\_input, complex\_t \* backward\_scratch, unsigned int flags=FFTW\_MEASURE) <br> |
|  int | [**nproc\_kspace**](classhacc_1_1Dfft.md#function-nproc-kspace-1-2) (int i) const<br> |
|  int const (& | [**nproc\_kspace**](classhacc_1_1Dfft.md#function-nproc-kspace-2-2) () <br> |
|  int | [**nproc\_rspace**](classhacc_1_1Dfft.md#function-nproc-rspace-1-2) (int i) const<br> |
|  int const (& | [**nproc\_rspace**](classhacc_1_1Dfft.md#function-nproc-rspace-2-2) () <br> |
|  MPI\_Comm | [**parent\_comm**](classhacc_1_1Dfft.md#function-parent-comm) () const<br> |
|  int | [**self\_kspace**](classhacc_1_1Dfft.md#function-self-kspace-1-2) (int i) const<br> |
|  int const (& | [**self\_kspace**](classhacc_1_1Dfft.md#function-self-kspace-2-2) () <br> |
|  int | [**self\_rspace**](classhacc_1_1Dfft.md#function-self-rspace-1-2) (int i) const<br> |
|  int const (& | [**self\_rspace**](classhacc_1_1Dfft.md#function-self-rspace-2-2) () <br> |
|   | [**~Dfft**](classhacc_1_1Dfft.md#function-dfft) () <br> |




## Protected Attributes

| Type | Name |
| ---: | :--- |
|  bool | [**PlansMade**](classhacc_1_1Dfft.md#variable-plansmade)  <br> |
|  Distribution & | [**d**](classhacc_1_1Dfft.md#variable-d)  <br> |
|  complex\_t \* | [**m\_bi**](classhacc_1_1Dfft.md#variable-m-bi)  <br> |
|  complex\_t \* | [**m\_bs**](classhacc_1_1Dfft.md#variable-m-bs)  <br> |
|  complex\_t \* | [**m\_fo**](classhacc_1_1Dfft.md#variable-m-fo)  <br> |
|  complex\_t \* | [**m\_fs**](classhacc_1_1Dfft.md#variable-m-fs)  <br> |
|  fftw\_plan | [**m\_plan\_b\_x**](classhacc_1_1Dfft.md#variable-m-plan-b-x)  <br> |
|  fftw\_plan | [**m\_plan\_b\_y**](classhacc_1_1Dfft.md#variable-m-plan-b-y)  <br> |
|  fftw\_plan | [**m\_plan\_b\_z**](classhacc_1_1Dfft.md#variable-m-plan-b-z)  <br> |
|  fftw\_plan | [**m\_plan\_f\_x**](classhacc_1_1Dfft.md#variable-m-plan-f-x)  <br> |
|  fftw\_plan | [**m\_plan\_f\_y**](classhacc_1_1Dfft.md#variable-m-plan-f-y)  <br> |
|  fftw\_plan | [**m\_plan\_f\_z**](classhacc_1_1Dfft.md#variable-m-plan-f-z)  <br> |




## Public Functions Documentation


### <a href="#function-dfft-1-2" id="function-dfft-1-2">function Dfft [1/2]</a>


```cpp
inline hacc::Dfft::Dfft (
    Distribution & dist
) 
```



### <a href="#function-dfft-2-2" id="function-dfft-2-2">function Dfft [2/2]</a>


```cpp
inline hacc::Dfft::Dfft (
    Distribution & dist,
    complex_t * forward_output,
    complex_t * forward_scratch,
    complex_t * backward_input,
    complex_t * backward_scratch,
    unsigned int flags=FFTW_MEASURE
) 
```



### <a href="#function-backward-1-3" id="function-backward-1-3">function backward [1/3]</a>


```cpp
inline void hacc::Dfft::backward (
    complex_t * out
) 
```



### <a href="#function-backward-2-3" id="function-backward-2-3">function backward [2/3]</a>


```cpp
inline void hacc::Dfft::backward (
    float * out,
    size_t ghost0,
    size_t ghost1
) 
```



### <a href="#function-backward-3-3" id="function-backward-3-3">function backward [3/3]</a>


```cpp
inline void hacc::Dfft::backward (
    float * out,
    size_t ghost
) 
```



### <a href="#function-cartcomm-kspace" id="function-cartcomm-kspace">function cartcomm\_kspace </a>


```cpp
inline MPI_Comm hacc::Dfft::cartcomm_kspace () const
```



### <a href="#function-cartcomm-rspace" id="function-cartcomm-rspace">function cartcomm\_rspace </a>


```cpp
inline MPI_Comm hacc::Dfft::cartcomm_rspace () const
```



### <a href="#function-forward-1-3" id="function-forward-1-3">function forward [1/3]</a>


```cpp
inline void hacc::Dfft::forward (
    complex_t const * in
) 
```



### <a href="#function-forward-2-3" id="function-forward-2-3">function forward [2/3]</a>


```cpp
inline void hacc::Dfft::forward (
    float const * in,
    size_t ghost0,
    size_t ghost1
) 
```



### <a href="#function-forward-3-3" id="function-forward-3-3">function forward [3/3]</a>


```cpp
inline void hacc::Dfft::forward (
    float const * in,
    size_t ghost
) 
```



### <a href="#function-get-d" id="function-get-d">function get\_d </a>


```cpp
inline Distribution & hacc::Dfft::get_d () 
```



### <a href="#function-global-ng-1-2" id="function-global-ng-1-2">function global\_ng [1/2]</a>


```cpp
inline int hacc::Dfft::global_ng (
    int i
) const
```



### <a href="#function-global-ng-2-2" id="function-global-ng-2-2">function global\_ng [2/2]</a>


```cpp
inline int const (& hacc::Dfft::global_ng () 
```



### <a href="#function-global-size" id="function-global-size">function global\_size </a>


```cpp
inline size_t hacc::Dfft::global_size () const
```



### <a href="#function-local-ng-kspace-1-2" id="function-local-ng-kspace-1-2">function local\_ng\_kspace [1/2]</a>


```cpp
inline int hacc::Dfft::local_ng_kspace (
    int i
) const
```



### <a href="#function-local-ng-kspace-2-2" id="function-local-ng-kspace-2-2">function local\_ng\_kspace [2/2]</a>


```cpp
inline int const (& hacc::Dfft::local_ng_kspace () 
```



### <a href="#function-local-ng-rspace-1-2" id="function-local-ng-rspace-1-2">function local\_ng\_rspace [1/2]</a>


```cpp
inline int hacc::Dfft::local_ng_rspace (
    int i
) const
```



### <a href="#function-local-ng-rspace-2-2" id="function-local-ng-rspace-2-2">function local\_ng\_rspace [2/2]</a>


```cpp
inline int const (& hacc::Dfft::local_ng_rspace () 
```



### <a href="#function-local-size" id="function-local-size">function local\_size </a>


```cpp
inline size_t hacc::Dfft::local_size () const
```



### <a href="#function-makeplans" id="function-makeplans">function makePlans </a>


```cpp
inline void hacc::Dfft::makePlans (
    complex_t * forward_output,
    complex_t * forward_scratch,
    complex_t * backward_input,
    complex_t * backward_scratch,
    unsigned int flags=FFTW_MEASURE
) 
```



### <a href="#function-nproc-kspace-1-2" id="function-nproc-kspace-1-2">function nproc\_kspace [1/2]</a>


```cpp
inline int hacc::Dfft::nproc_kspace (
    int i
) const
```



### <a href="#function-nproc-kspace-2-2" id="function-nproc-kspace-2-2">function nproc\_kspace [2/2]</a>


```cpp
inline int const (& hacc::Dfft::nproc_kspace () 
```



### <a href="#function-nproc-rspace-1-2" id="function-nproc-rspace-1-2">function nproc\_rspace [1/2]</a>


```cpp
inline int hacc::Dfft::nproc_rspace (
    int i
) const
```



### <a href="#function-nproc-rspace-2-2" id="function-nproc-rspace-2-2">function nproc\_rspace [2/2]</a>


```cpp
inline int const (& hacc::Dfft::nproc_rspace () 
```



### <a href="#function-parent-comm" id="function-parent-comm">function parent\_comm </a>


```cpp
inline MPI_Comm hacc::Dfft::parent_comm () const
```



### <a href="#function-self-kspace-1-2" id="function-self-kspace-1-2">function self\_kspace [1/2]</a>


```cpp
inline int hacc::Dfft::self_kspace (
    int i
) const
```



### <a href="#function-self-kspace-2-2" id="function-self-kspace-2-2">function self\_kspace [2/2]</a>


```cpp
inline int const (& hacc::Dfft::self_kspace () 
```



### <a href="#function-self-rspace-1-2" id="function-self-rspace-1-2">function self\_rspace [1/2]</a>


```cpp
inline int hacc::Dfft::self_rspace (
    int i
) const
```



### <a href="#function-self-rspace-2-2" id="function-self-rspace-2-2">function self\_rspace [2/2]</a>


```cpp
inline int const (& hacc::Dfft::self_rspace () 
```



### <a href="#function-dfft" id="function-dfft">function ~Dfft </a>


```cpp
inline hacc::Dfft::~Dfft () 
```


## Protected Attributes Documentation


### <a href="#variable-plansmade" id="variable-plansmade">variable PlansMade </a>


```cpp
bool hacc::Dfft::PlansMade;
```



### <a href="#variable-d" id="variable-d">variable d </a>


```cpp
Distribution& hacc::Dfft::d;
```



### <a href="#variable-m-bi" id="variable-m-bi">variable m\_bi </a>


```cpp
complex_t* hacc::Dfft::m_bi;
```



### <a href="#variable-m-bs" id="variable-m-bs">variable m\_bs </a>


```cpp
complex_t* hacc::Dfft::m_bs;
```



### <a href="#variable-m-fo" id="variable-m-fo">variable m\_fo </a>


```cpp
complex_t* hacc::Dfft::m_fo;
```



### <a href="#variable-m-fs" id="variable-m-fs">variable m\_fs </a>


```cpp
complex_t* hacc::Dfft::m_fs;
```



### <a href="#variable-m-plan-b-x" id="variable-m-plan-b-x">variable m\_plan\_b\_x </a>


```cpp
fftw_plan hacc::Dfft::m_plan_b_x;
```



### <a href="#variable-m-plan-b-y" id="variable-m-plan-b-y">variable m\_plan\_b\_y </a>


```cpp
fftw_plan hacc::Dfft::m_plan_b_y;
```



### <a href="#variable-m-plan-b-z" id="variable-m-plan-b-z">variable m\_plan\_b\_z </a>


```cpp
fftw_plan hacc::Dfft::m_plan_b_z;
```



### <a href="#variable-m-plan-f-x" id="variable-m-plan-f-x">variable m\_plan\_f\_x </a>


```cpp
fftw_plan hacc::Dfft::m_plan_f_x;
```



### <a href="#variable-m-plan-f-y" id="variable-m-plan-f-y">variable m\_plan\_f\_y </a>


```cpp
fftw_plan hacc::Dfft::m_plan_f_y;
```



### <a href="#variable-m-plan-f-z" id="variable-m-plan-f-z">variable m\_plan\_f\_z </a>


```cpp
fftw_plan hacc::Dfft::m_plan_f_z;
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Exec/Test_Only_Axions/Dfft.H`