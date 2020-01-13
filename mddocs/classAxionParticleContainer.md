
# Class AxionParticleContainer

**template &lt;int N&gt;**


[**Class List**](annotated.md) **>** [**AxionParticleContainer**](classAxionParticleContainer.md)








Inherits the following classes: ParticleContainer< N >








## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::deque&lt; ParticleType &gt; | [**PBox**](classAxionParticleContainer.md#typedef-pbox)  <br> |
| typedef std::map&lt; int, PBox &gt; | [**PMap**](classAxionParticleContainer.md#typedef-pmap)  <br> |
| typedef Particle&lt; N &gt; | [**ParticleType**](classAxionParticleContainer.md#typedef-particletype)  <br> |



## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  Real | [**damping\_constant**](classAxionParticleContainer.md#variable-damping-constant)   = = 0.0<br> |
|  bool | [**initzd**](classAxionParticleContainer.md#variable-initzd)   = = false<br> |
|  Real | [**max\_acceleration**](classAxionParticleContainer.md#variable-max-acceleration)   = = 1e20<br> |
|  MultiFab \* | [**mf\_pointer**](classAxionParticleContainer.md#variable-mf-pointer)   = = 0<br> |
|  bool | [**particles\_moved**](classAxionParticleContainer.md#variable-particles-moved)   = = true<br> |
|  int | [**smoothing\_length**](classAxionParticleContainer.md#variable-smoothing-length)   = = 4<br> |

## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**AssignDensity**](classAxionParticleContainer.md#function-assigndensity) (PArray&lt; MultiFab &gt; & mf, int lev\_min=0, int ncomp=1, int finest\_level=-1, double masstreshold=1.e50) const<br> |
|  void | [**AssignDensityGradientsSingleLevel**](classAxionParticleContainer.md#function-assigndensitygradientssinglelevel) (MultiFab & mf\_to\_be\_filled, int lev) const<br> |
|  void | [**AssignDensitySingleLevel**](classAxionParticleContainer.md#function-assigndensitysinglelevel) (MultiFab & mf, int level, int ncomp=1, int particle\_lvl\_offset=0) const<br> |
|   | [**AxionParticleContainer**](classAxionParticleContainer.md#function-axionparticlecontainer) (Amr \* amr) <br> |
|  void | [**InitEmpty**](classAxionParticleContainer.md#function-initempty) (const int max\_level) <br> |
|  void | [**InitVarCount**](classAxionParticleContainer.md#function-initvarcount) (MultiFab & mf, long n\_axpart) <br> |
|  void | [**InitVarMass**](classAxionParticleContainer.md#function-initvarmass) (MultiFab & mf) <br> |
|  void | [**addArtificialViscosity**](classAxionParticleContainer.md#function-addartificialviscosity) (MultiFab & acc\_vector, int lev, Real a\_old, Real a\_half) <br> |
|  void | [**addQuantumPressure**](classAxionParticleContainer.md#function-addquantumpressure) (MultiFab & acc\_vector, int lev, Real a\_old, Real a\_half) <br> |
|  void | [**averageVelocity**](classAxionParticleContainer.md#function-averagevelocity) (int level) <br> |
|  void | [**computeEnergies**](classAxionParticleContainer.md#function-computeenergies) (Real cur\_time, MultiFab & Phi, Real a) <br> |
|  void | [**computeQuantumPressure**](classAxionParticleContainer.md#function-computequantumpressure) (MultiFab & axPressForceMf, MultiFab & pressureMf, int lev) <br> |
|  void | [**doFiltering**](classAxionParticleContainer.md#function-dofiltering) (MultiFab & mf, int comp, int niter, int nstride) const<br> |
|  Real | [**estTimestep**](classAxionParticleContainer.md#function-esttimestep) (MultiFab & grav\_vector, Real a, int lev, Real cfl) <br> |
|  void | [**moveKick**](classAxionParticleContainer.md#function-movekick) (MultiFab & grav\_vector, int level, Real timestep, Real a\_new=1.0, Real a\_half=1.0, int start\_comp\_for\_accel=-1) <br> |
|  void | [**moveKickDrift**](classAxionParticleContainer.md#function-movekickdrift) (MultiFab & grav\_vector, int level, Real timestep, Real a\_old=1.0, Real a\_half=1.0) <br> |
|  void | [**quantummoveKick**](classAxionParticleContainer.md#function-quantummovekick) (Real timestep, Real a\_old, Real a\_half, bool before\_drift) <br> |
|  void | [**remap**](classAxionParticleContainer.md#function-remap) () <br> |








## Public Types Documentation


### <a href="#typedef-pbox" id="typedef-pbox">typedef PBox </a>


```cpp
typedef std::deque<ParticleType> AxionParticleContainer< N >::PBox;
```



### <a href="#typedef-pmap" id="typedef-pmap">typedef PMap </a>


```cpp
typedef std::map<int,PBox> AxionParticleContainer< N >::PMap;
```



### <a href="#typedef-particletype" id="typedef-particletype">typedef ParticleType </a>


```cpp
typedef Particle<N> AxionParticleContainer< N >::ParticleType;
```


## Public Static Attributes Documentation


### <a href="#variable-damping-constant" id="variable-damping-constant">variable damping\_constant </a>


```cpp
Real AxionParticleContainer< N >::damping_constant;
```



### <a href="#variable-initzd" id="variable-initzd">variable initzd </a>


```cpp
bool AxionParticleContainer< N >::initzd;
```



### <a href="#variable-max-acceleration" id="variable-max-acceleration">variable max\_acceleration </a>


```cpp
Real AxionParticleContainer< N >::max_acceleration;
```



### <a href="#variable-mf-pointer" id="variable-mf-pointer">variable mf\_pointer </a>


```cpp
MultiFab * AxionParticleContainer< N >::mf_pointer;
```



### <a href="#variable-particles-moved" id="variable-particles-moved">variable particles\_moved </a>


```cpp
bool AxionParticleContainer< N >::particles_moved;
```



### <a href="#variable-smoothing-length" id="variable-smoothing-length">variable smoothing\_length </a>


```cpp
int AxionParticleContainer< N >::smoothing_length;
```


## Public Functions Documentation


### <a href="#function-assigndensity" id="function-assigndensity">function AssignDensity </a>


```cpp
void AxionParticleContainer::AssignDensity (
    PArray< MultiFab > & mf,
    int lev_min=0,
    int ncomp=1,
    int finest_level=-1,
    double masstreshold=1.e50
) const
```



### <a href="#function-assigndensitygradientssinglelevel" id="function-assigndensitygradientssinglelevel">function AssignDensityGradientsSingleLevel </a>


```cpp
void AxionParticleContainer::AssignDensityGradientsSingleLevel (
    MultiFab & mf_to_be_filled,
    int lev
) const
```



### <a href="#function-assigndensitysinglelevel" id="function-assigndensitysinglelevel">function AssignDensitySingleLevel </a>


```cpp
void AxionParticleContainer::AssignDensitySingleLevel (
    MultiFab & mf,
    int level,
    int ncomp=1,
    int particle_lvl_offset=0
) const
```



### <a href="#function-axionparticlecontainer" id="function-axionparticlecontainer">function AxionParticleContainer </a>


```cpp
inline AxionParticleContainer::AxionParticleContainer (
    Amr * amr
) 
```



### <a href="#function-initempty" id="function-initempty">function InitEmpty </a>


```cpp
void AxionParticleContainer::InitEmpty (
    const int max_level
) 
```



### <a href="#function-initvarcount" id="function-initvarcount">function InitVarCount </a>


```cpp
void AxionParticleContainer::InitVarCount (
    MultiFab & mf,
    long n_axpart
) 
```



### <a href="#function-initvarmass" id="function-initvarmass">function InitVarMass </a>


```cpp
void AxionParticleContainer::InitVarMass (
    MultiFab & mf
) 
```



### <a href="#function-addartificialviscosity" id="function-addartificialviscosity">function addArtificialViscosity </a>


```cpp
void AxionParticleContainer::addArtificialViscosity (
    MultiFab & acc_vector,
    int lev,
    Real a_old,
    Real a_half
) 
```



### <a href="#function-addquantumpressure" id="function-addquantumpressure">function addQuantumPressure </a>


```cpp
void AxionParticleContainer::addQuantumPressure (
    MultiFab & acc_vector,
    int lev,
    Real a_old,
    Real a_half
) 
```



### <a href="#function-averagevelocity" id="function-averagevelocity">function averageVelocity </a>


```cpp
void AxionParticleContainer::averageVelocity (
    int level
) 
```



### <a href="#function-computeenergies" id="function-computeenergies">function computeEnergies </a>


```cpp
void AxionParticleContainer::computeEnergies (
    Real cur_time,
    MultiFab & Phi,
    Real a
) 
```



### <a href="#function-computequantumpressure" id="function-computequantumpressure">function computeQuantumPressure </a>


```cpp
void AxionParticleContainer::computeQuantumPressure (
    MultiFab & axPressForceMf,
    MultiFab & pressureMf,
    int lev
) 
```



### <a href="#function-dofiltering" id="function-dofiltering">function doFiltering </a>


```cpp
void AxionParticleContainer::doFiltering (
    MultiFab & mf,
    int comp,
    int niter,
    int nstride
) const
```



### <a href="#function-esttimestep" id="function-esttimestep">function estTimestep </a>


```cpp
Real AxionParticleContainer::estTimestep (
    MultiFab & grav_vector,
    Real a,
    int lev,
    Real cfl
) 
```



### <a href="#function-movekick" id="function-movekick">function moveKick </a>


```cpp
void AxionParticleContainer::moveKick (
    MultiFab & grav_vector,
    int level,
    Real timestep,
    Real a_new=1.0,
    Real a_half=1.0,
    int start_comp_for_accel=-1
) 
```



### <a href="#function-movekickdrift" id="function-movekickdrift">function moveKickDrift </a>


```cpp
void AxionParticleContainer::moveKickDrift (
    MultiFab & grav_vector,
    int level,
    Real timestep,
    Real a_old=1.0,
    Real a_half=1.0
) 
```



### <a href="#function-quantummovekick" id="function-quantummovekick">function quantummoveKick </a>


```cpp
void AxionParticleContainer::quantummoveKick (
    Real timestep,
    Real a_old,
    Real a_half,
    bool before_drift
) 
```



### <a href="#function-remap" id="function-remap">function remap </a>


```cpp
void AxionParticleContainer::remap () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/AxParticles.H`