
# Class StochasticForcing


[**Class List**](annotated.md) **>** [**StochasticForcing**](classStochasticForcing.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**ReadSpectrum**](classStochasticForcing.md#function-readspectrum) (char \* fname) <br> |
|   | [**StochasticForcing**](classStochasticForcing.md#function-stochasticforcing) () <br> |
|  int | [**WriteSpectrum**](classStochasticForcing.md#function-writespectrum) (char \* fname) <br> |
|  void | [**copy\_ExpandedSpectrum**](classStochasticForcing.md#function-copy-expandedspectrum) (int dim, amrex::Real \* target) <br> |
|  void | [**copy\_SpectrumEven**](classStochasticForcing.md#function-copy-spectrumeven) (int dim, amrex::Real \* target) <br> |
|  void | [**copy\_SpectrumOdd**](classStochasticForcing.md#function-copy-spectrumodd) (int dim, amrex::Real \* target) <br> |
|  void | [**copy\_mask**](classStochasticForcing.md#function-copy-mask) (int \* target) <br> |
|  void | [**distribute**](classStochasticForcing.md#function-distribute) (void) <br> |
|  void | [**evolve**](classStochasticForcing.md#function-evolve) (amrex::Real dt) <br> |
|  amrex::Real | [**get\_IntgrLength**](classStochasticForcing.md#function-get-intgrlength) (int dim) <br> |
|  amrex::Real | [**get\_IntgrTime**](classStochasticForcing.md#function-get-intgrtime) (int dim) <br> |
|  int | [**get\_LeftBoundary**](classStochasticForcing.md#function-get-leftboundary) (int dim) <br> |
|  int | [**get\_NumModes**](classStochasticForcing.md#function-get-nummodes) (void) <br> |
|  int | [**get\_NumNonZeroModes**](classStochasticForcing.md#function-get-numnonzeromodes) (void) <br> |
|  int | [**get\_Range**](classStochasticForcing.md#function-get-range) (int dim) <br> |
|  int | [**get\_RightBoundary**](classStochasticForcing.md#function-get-rightboundary) (int dim) <br> |
|  spect\_profile\_type | [**get\_SpectProfile**](classStochasticForcing.md#function-get-spectprofile) (void) <br> |
|  int | [**get\_SpectralRank**](classStochasticForcing.md#function-get-spectralrank) (void) <br> |
|  amrex::Real | [**get\_WaveNumber**](classStochasticForcing.md#function-get-wavenumber) (int dim) <br> |
|  void | [**init**](classStochasticForcing.md#function-init) (int rank, const amrex::Real \* prob\_lo, const amrex::Real \* prob\_hi) <br> |
|  void | [**read\_Spectrum**](classStochasticForcing.md#function-read-spectrum) (std::ifstream & input) <br> |
|  amrex::Real | [**rms**](classStochasticForcing.md#function-rms) (void) <br> |
|  void | [**set\_SolenoidalWeight**](classStochasticForcing.md#function-set-solenoidalweight) (int my\_soln\_weight) <br> |
|  void | [**set\_decay**](classStochasticForcing.md#function-set-decay) (void) <br> |
|  void | [**write\_Spectrum**](classStochasticForcing.md#function-write-spectrum) (std::ofstream & output) <br> |
|   | [**~StochasticForcing**](classStochasticForcing.md#function-stochasticforcing) () <br> |








## Public Functions Documentation


### <a href="#function-readspectrum" id="function-readspectrum">function ReadSpectrum </a>


```cpp
int StochasticForcing::ReadSpectrum (
    char * fname
) 
```



### <a href="#function-stochasticforcing" id="function-stochasticforcing">function StochasticForcing </a>


```cpp
StochasticForcing::StochasticForcing () 
```



### <a href="#function-writespectrum" id="function-writespectrum">function WriteSpectrum </a>


```cpp
int StochasticForcing::WriteSpectrum (
    char * fname
) 
```



### <a href="#function-copy-expandedspectrum" id="function-copy-expandedspectrum">function copy\_ExpandedSpectrum </a>


```cpp
inline void StochasticForcing::copy_ExpandedSpectrum (
    int dim,
    amrex::Real * target
) 
```



### <a href="#function-copy-spectrumeven" id="function-copy-spectrumeven">function copy\_SpectrumEven </a>


```cpp
inline void StochasticForcing::copy_SpectrumEven (
    int dim,
    amrex::Real * target
) 
```



### <a href="#function-copy-spectrumodd" id="function-copy-spectrumodd">function copy\_SpectrumOdd </a>


```cpp
inline void StochasticForcing::copy_SpectrumOdd (
    int dim,
    amrex::Real * target
) 
```



### <a href="#function-copy-mask" id="function-copy-mask">function copy\_mask </a>


```cpp
inline void StochasticForcing::copy_mask (
    int * target
) 
```



### <a href="#function-distribute" id="function-distribute">function distribute </a>


```cpp
void StochasticForcing::distribute (
    void
) 
```



### <a href="#function-evolve" id="function-evolve">function evolve </a>


```cpp
void StochasticForcing::evolve (
    amrex::Real dt
) 
```



### <a href="#function-get-intgrlength" id="function-get-intgrlength">function get\_IntgrLength </a>


```cpp
inline amrex::Real StochasticForcing::get_IntgrLength (
    int dim
) 
```



### <a href="#function-get-intgrtime" id="function-get-intgrtime">function get\_IntgrTime </a>


```cpp
inline amrex::Real StochasticForcing::get_IntgrTime (
    int dim
) 
```



### <a href="#function-get-leftboundary" id="function-get-leftboundary">function get\_LeftBoundary </a>


```cpp
inline int StochasticForcing::get_LeftBoundary (
    int dim
) 
```



### <a href="#function-get-nummodes" id="function-get-nummodes">function get\_NumModes </a>


```cpp
inline int StochasticForcing::get_NumModes (
    void
) 
```



### <a href="#function-get-numnonzeromodes" id="function-get-numnonzeromodes">function get\_NumNonZeroModes </a>


```cpp
inline int StochasticForcing::get_NumNonZeroModes (
    void
) 
```



### <a href="#function-get-range" id="function-get-range">function get\_Range </a>


```cpp
int StochasticForcing::get_Range (
    int dim
) 
```



### <a href="#function-get-rightboundary" id="function-get-rightboundary">function get\_RightBoundary </a>


```cpp
inline int StochasticForcing::get_RightBoundary (
    int dim
) 
```



### <a href="#function-get-spectprofile" id="function-get-spectprofile">function get\_SpectProfile </a>


```cpp
inline spect_profile_type StochasticForcing::get_SpectProfile (
    void
) 
```



### <a href="#function-get-spectralrank" id="function-get-spectralrank">function get\_SpectralRank </a>


```cpp
inline int StochasticForcing::get_SpectralRank (
    void
) 
```



### <a href="#function-get-wavenumber" id="function-get-wavenumber">function get\_WaveNumber </a>


```cpp
inline amrex::Real StochasticForcing::get_WaveNumber (
    int dim
) 
```



### <a href="#function-init" id="function-init">function init </a>


```cpp
void StochasticForcing::init (
    int rank,
    const amrex::Real * prob_lo,
    const amrex::Real * prob_hi
) 
```



### <a href="#function-read-spectrum" id="function-read-spectrum">function read\_Spectrum </a>


```cpp
inline void StochasticForcing::read_Spectrum (
    std::ifstream & input
) 
```



### <a href="#function-rms" id="function-rms">function rms </a>


```cpp
amrex::Real StochasticForcing::rms (
    void
) 
```



### <a href="#function-set-solenoidalweight" id="function-set-solenoidalweight">function set\_SolenoidalWeight </a>


```cpp
inline void StochasticForcing::set_SolenoidalWeight (
    int my_soln_weight
) 
```



### <a href="#function-set-decay" id="function-set-decay">function set\_decay </a>


```cpp
inline void StochasticForcing::set_decay (
    void
) 
```



### <a href="#function-write-spectrum" id="function-write-spectrum">function write\_Spectrum </a>


```cpp
inline void StochasticForcing::write_Spectrum (
    std::ofstream & output
) 
```



### <a href="#function-stochasticforcing" id="function-stochasticforcing">function ~StochasticForcing </a>


```cpp
StochasticForcing::~StochasticForcing () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/Forcing/Forcing.H`