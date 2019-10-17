
# Class NyxParticleContainerBase


[**Class List**](annotated.md) **>** [**NyxParticleContainerBase**](classNyxParticleContainerBase.md)










Inherited by the following classes: [NyxParticleContainer](classNyxParticleContainer.md),  [NyxParticleContainer](classNyxParticleContainer.md),  [NyxParticleContainer](classNyxParticleContainer.md)










## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**AssignDensity**](classNyxParticleContainerBase.md#function-assigndensity) (amrex::Vector&lt; std::unique\_ptr&lt; amrex::MultiFab &gt; &gt; & mf, int lev\_min=0, int ncomp=1, int finest\_level=-1, int ngrow=1) const = 0<br> |
| virtual void | [**AssignDensitySingleLevel**](classNyxParticleContainerBase.md#function-assigndensitysinglelevel) (amrex::MultiFab & mf, int level, int ncomp=1, int particle\_lvl\_offset=0) const = 0<br> |
| virtual void | [**Redistribute**](classNyxParticleContainerBase.md#function-redistribute) (int lev\_min=0, int lev\_max=-1, int nGrow=0) = 0<br> |
| virtual void | [**RemoveParticlesAtLevel**](classNyxParticleContainerBase.md#function-removeparticlesatlevel) (int level) = 0<br> |
| virtual int | [**finestLevel**](classNyxParticleContainerBase.md#function-finestlevel) () const = 0<br> |
| virtual void | [**moveKick**](classNyxParticleContainerBase.md#function-movekick) (amrex::MultiFab & acceleration, int level, amrex::Real timestep, amrex::Real a\_new=1.0, amrex::Real a\_half=1.0) = 0<br> |
| virtual void | [**moveKickDrift**](classNyxParticleContainerBase.md#function-movekickdrift) (amrex::MultiFab & acceleration, int level, amrex::Real timestep, amrex::Real a\_old=1.0, amrex::Real a\_half=1.0, int where\_width=0) = 0<br> |
| virtual amrex::Real | [**sumParticleMass**](classNyxParticleContainerBase.md#function-sumparticlemass) (int level) const = 0<br> |
| virtual  | [**~NyxParticleContainerBase**](classNyxParticleContainerBase.md#function-nyxparticlecontainerbase) () <br> |








## Public Functions Documentation


### <a href="#function-assigndensity" id="function-assigndensity">function AssignDensity </a>


```cpp
virtual void NyxParticleContainerBase::AssignDensity (
    amrex::Vector< std::unique_ptr< amrex::MultiFab > > & mf,
    int lev_min=0,
    int ncomp=1,
    int finest_level=-1,
    int ngrow=1
) const = 0
```



### <a href="#function-assigndensitysinglelevel" id="function-assigndensitysinglelevel">function AssignDensitySingleLevel </a>


```cpp
virtual void NyxParticleContainerBase::AssignDensitySingleLevel (
    amrex::MultiFab & mf,
    int level,
    int ncomp=1,
    int particle_lvl_offset=0
) const = 0
```



### <a href="#function-redistribute" id="function-redistribute">function Redistribute </a>


```cpp
virtual void NyxParticleContainerBase::Redistribute (
    int lev_min=0,
    int lev_max=-1,
    int nGrow=0
) = 0
```



### <a href="#function-removeparticlesatlevel" id="function-removeparticlesatlevel">function RemoveParticlesAtLevel </a>


```cpp
virtual void NyxParticleContainerBase::RemoveParticlesAtLevel (
    int level
) = 0
```



### <a href="#function-finestlevel" id="function-finestlevel">function finestLevel </a>


```cpp
virtual int NyxParticleContainerBase::finestLevel () const = 0
```



### <a href="#function-movekick" id="function-movekick">function moveKick </a>


```cpp
virtual void NyxParticleContainerBase::moveKick (
    amrex::MultiFab & acceleration,
    int level,
    amrex::Real timestep,
    amrex::Real a_new=1.0,
    amrex::Real a_half=1.0
) = 0
```



### <a href="#function-movekickdrift" id="function-movekickdrift">function moveKickDrift </a>


```cpp
virtual void NyxParticleContainerBase::moveKickDrift (
    amrex::MultiFab & acceleration,
    int level,
    amrex::Real timestep,
    amrex::Real a_old=1.0,
    amrex::Real a_half=1.0,
    int where_width=0
) = 0
```



### <a href="#function-sumparticlemass" id="function-sumparticlemass">function sumParticleMass </a>


```cpp
virtual amrex::Real NyxParticleContainerBase::sumParticleMass (
    int level
) const = 0
```



### <a href="#function-nyxparticlecontainerbase" id="function-nyxparticlecontainerbase">function ~NyxParticleContainerBase </a>


```cpp
inline virtual NyxParticleContainerBase::~NyxParticleContainerBase () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/NyxParticleContainer.H`