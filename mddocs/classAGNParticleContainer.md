
# Class AGNParticleContainer


[**Class List**](annotated.md) **>** [**AGNParticleContainer**](classAGNParticleContainer.md)








Inherits the following classes: [NyxParticleContainer](classNyxParticleContainer.md)








## Public Types

| Type | Name |
| ---: | :--- |
| typedef amrex::ParIter&lt; 3+BL\_SPACEDIM &gt; | [**MyParIter**](classAGNParticleContainer.md#typedef-mypariter)  <br> |

## Public Types inherited from [NyxParticleContainer](classNyxParticleContainer.md)

| Type | Name |
| ---: | :--- |
| typedef typename amrex::ParticleContainer&lt; NSR, NSI, NAR, NAI &gt;::AoS | [**AoS**](classNyxParticleContainer.md#typedef-aos)  <br> |
| typedef amrex::ParConstIter&lt; NSR, NSI, NAR, NAI &gt; | [**MyConstParIter**](classNyxParticleContainer.md#typedef-myconstpariter)  <br> |
| typedef amrex::ParIter&lt; NSR, NSI, NAR, NAI &gt; | [**MyParIter**](classNyxParticleContainer.md#typedef-mypariter)  <br> |
| typedef typename amrex::ParticleContainer&lt; NSR, NSI, NAR, NAI &gt;::ParticleLevel | [**ParticleLevel**](classNyxParticleContainer.md#typedef-particlelevel)  <br> |
| typedef amrex::ParticleTile&lt; NSR, NSI, NAR, NAI &gt; | [**ParticleTileType**](classNyxParticleContainer.md#typedef-particletiletype)  <br> |
| typedef amrex::Particle&lt; NSR, NSI &gt; | [**ParticleType**](classNyxParticleContainer.md#typedef-particletype)  <br> |











## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**AGNParticleContainer**](classAGNParticleContainer.md#function-agnparticlecontainer) (amrex::Amr \* amr, int nghost) <br> |
|  void | [**AccreteMass**](classAGNParticleContainer.md#function-accretemass) (int lev, amrex::MultiFab & state, amrex::MultiFab & density\_lost, amrex::Real dt) <br> |
|  void | [**AddOneParticle**](classAGNParticleContainer.md#function-addoneparticle-1-2) (int lev, int grid, int tile, amrex::Real mass, amrex::Real x, amrex::Real y, amrex::Real z) <br> |
|  void | [**AddOneParticle**](classAGNParticleContainer.md#function-addoneparticle-2-2) (ParticleTileType & particle\_tile, amrex::Real mass, amrex::Real x, amrex::Real y, amrex::Real z) <br> |
|  void | [**ComputeOverlap**](classAGNParticleContainer.md#function-computeoverlap) (int lev) <br> |
|  void | [**ComputeParticleVelocity**](classAGNParticleContainer.md#function-computeparticlevelocity) (int lev, amrex::MultiFab & state\_old, amrex::MultiFab & state\_new, int add\_energy) <br> |
|  void | [**Merge**](classAGNParticleContainer.md#function-merge) (int lev) <br> |
|  const int | [**NumberOfParticles**](classAGNParticleContainer.md#function-numberofparticles) (MyParIter & pti) <br> |
|  void | [**ReleaseEnergy**](classAGNParticleContainer.md#function-releaseenergy) (int lev, amrex::MultiFab & state, amrex::MultiFab & D\_new, amrex::Real a) <br> |
| virtual void | [**moveKick**](classAGNParticleContainer.md#function-movekick) (amrex::MultiFab & acceleration, int level, amrex::Real timestep, amrex::Real a\_new=1.0, amrex::Real a\_half=1.0) <br> |
| virtual void | [**moveKickDrift**](classAGNParticleContainer.md#function-movekickdrift) (amrex::MultiFab & acceleration, int level, amrex::Real timestep, amrex::Real a\_old=1.0, amrex::Real a\_half=1.0, int where\_width=0) <br> |
|  void | [**writeAllAtLevel**](classAGNParticleContainer.md#function-writeallatlevel) (int lev) <br> |
| virtual  | [**~AGNParticleContainer**](classAGNParticleContainer.md#function-agnparticlecontainer) () <br> |

## Public Functions inherited from [NyxParticleContainer](classNyxParticleContainer.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**AssignDensity**](classNyxParticleContainer.md#function-assigndensity) (amrex::Vector&lt; std::unique\_ptr&lt; amrex::MultiFab &gt; &gt; & mf, int lev\_min=0, int ncomp=1, int finest\_level=-1, int ngrow=1) override const<br> |
| virtual void | [**AssignDensitySingleLevel**](classNyxParticleContainer.md#function-assigndensitysinglelevel) (amrex::MultiFab & mf, int level, int ncomp=1, int particle\_lvl\_offset=0) override const<br> |
|  void | [**GetParticleVelocities**](classNyxParticleContainer.md#function-getparticlevelocities) (amrex::Vector&lt; amrex::Real &gt; & part\_vels) <br> |
|  void | [**MultiplyParticleMass**](classNyxParticleContainer.md#function-multiplyparticlemass) (int lev, amrex::Real mult) <br> |
| virtual void | [**NyxCheckpoint**](classNyxParticleContainer.md#function-nyxcheckpoint) (const std::string & dir, const std::string & name) const<br> |
|   | [**NyxParticleContainer**](classNyxParticleContainer.md#function-nyxparticlecontainer) (amrex::Amr \* amr, int nghost=0) <br> |
| virtual void | [**Redistribute**](classNyxParticleContainer.md#function-redistribute) (int lev\_min=0, int lev\_max=-1, int nGrow=0) override<br> |
| virtual void | [**RemoveParticlesAtLevel**](classNyxParticleContainer.md#function-removeparticlesatlevel) (int level) override<br> |
|  void | [**SetParticleVelocities**](classNyxParticleContainer.md#function-setparticlevelocities) (amrex::Vector&lt; amrex::Real &gt; & part\_data) <br> |
| virtual void | [**WriteNyxPlotFile**](classNyxParticleContainer.md#function-writenyxplotfile) (const std::string & dir, const std::string & name) const<br> |
|  amrex::Real | [**estTimestep**](classNyxParticleContainer.md#function-esttimestep-1-2) (amrex::MultiFab & acceleration, int level, amrex::Real cfl) const<br> |
|  amrex::Real | [**estTimestep**](classNyxParticleContainer.md#function-esttimestep-2-2) (amrex::MultiFab & acceleration, amrex::Real a, int level, amrex::Real cfl) const<br> |
| virtual int | [**finestLevel**](classNyxParticleContainer.md#function-finestlevel) () override const<br> |
| virtual amrex::Real | [**sumParticleMass**](classNyxParticleContainer.md#function-sumparticlemass) (int level) override const<br> |
|  void | [**sumParticleMomentum**](classNyxParticleContainer.md#function-sumparticlemomentum) (int lev, amrex::Real \* mom) const<br> |
| virtual  | [**~NyxParticleContainer**](classNyxParticleContainer.md#function-nyxparticlecontainer) () <br> |

## Public Functions inherited from [NyxParticleContainerBase](classNyxParticleContainerBase.md)

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










## Protected Attributes

| Type | Name |
| ---: | :--- |
|  amrex::Vector&lt; std::string &gt; | [**real\_comp\_names**](classAGNParticleContainer.md#variable-real-comp-names)  <br> |
|  bool | [**sub\_cycle**](classAGNParticleContainer.md#variable-sub-cycle)  <br> |

## Protected Attributes inherited from [NyxParticleContainer](classNyxParticleContainer.md)

| Type | Name |
| ---: | :--- |
|  amrex::Vector&lt; std::string &gt; | [**real\_comp\_names**](classNyxParticleContainer.md#variable-real-comp-names)  <br> |
|  bool | [**sub\_cycle**](classNyxParticleContainer.md#variable-sub-cycle)  <br> |











## Public Types Documentation


### <a href="#typedef-mypariter" id="typedef-mypariter">typedef MyParIter </a>


```cpp
using AGNParticleContainer::MyParIter =  amrex::ParIter<3+BL_SPACEDIM>;
```


## Public Functions Documentation


### <a href="#function-agnparticlecontainer" id="function-agnparticlecontainer">function AGNParticleContainer </a>


```cpp
inline AGNParticleContainer::AGNParticleContainer (
    amrex::Amr * amr,
    int nghost
) 
```



### <a href="#function-accretemass" id="function-accretemass">function AccreteMass </a>


```cpp
void AGNParticleContainer::AccreteMass (
    int lev,
    amrex::MultiFab & state,
    amrex::MultiFab & density_lost,
    amrex::Real dt
) 
```


Accrete mass from the grid onto the existing AGN particles 


        

### <a href="#function-addoneparticle-1-2" id="function-addoneparticle-1-2">function AddOneParticle [1/2]</a>


```cpp
inline void AGNParticleContainer::AddOneParticle (
    int lev,
    int grid,
    int tile,
    amrex::Real mass,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z
) 
```



### <a href="#function-addoneparticle-2-2" id="function-addoneparticle-2-2">function AddOneParticle [2/2]</a>


```cpp
inline void AGNParticleContainer::AddOneParticle (
    ParticleTileType & particle_tile,
    amrex::Real mass,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z
) 
```



### <a href="#function-computeoverlap" id="function-computeoverlap">function ComputeOverlap </a>


```cpp
void AGNParticleContainer::ComputeOverlap (
    int lev
) 
```


Invalidate particles in cells that are already occupied 


        

### <a href="#function-computeparticlevelocity" id="function-computeparticlevelocity">function ComputeParticleVelocity </a>


```cpp
void AGNParticleContainer::ComputeParticleVelocity (
    int lev,
    amrex::MultiFab & state_old,
    amrex::MultiFab & state_new,
    int add_energy
) 
```


Compute the momentum that has been removed from the gas in order to define the particle velocity 


        

### <a href="#function-merge" id="function-merge">function Merge </a>


```cpp
void AGNParticleContainer::Merge (
    int lev
) 
```


Invalidate particles that have been merged with other particles 


        

### <a href="#function-numberofparticles" id="function-numberofparticles">function NumberOfParticles </a>


```cpp
inline const int AGNParticleContainer::NumberOfParticles (
    MyParIter & pti
) 
```



### <a href="#function-releaseenergy" id="function-releaseenergy">function ReleaseEnergy </a>


```cpp
void AGNParticleContainer::ReleaseEnergy (
    int lev,
    amrex::MultiFab & state,
    amrex::MultiFab & D_new,
    amrex::Real a
) 
```


Release energy if it exceeds thermal feedback threshold. 


        

### <a href="#function-movekick" id="function-movekick">function moveKick </a>


```cpp
virtual void AGNParticleContainer::moveKick (
    amrex::MultiFab & acceleration,
    int level,
    amrex::Real timestep,
    amrex::Real a_new=1.0,
    amrex::Real a_half=1.0
) 
```


Implements [*NyxParticleContainerBase::moveKick*](classNyxParticleContainerBase.md#function-movekick)


### <a href="#function-movekickdrift" id="function-movekickdrift">function moveKickDrift </a>


```cpp
virtual void AGNParticleContainer::moveKickDrift (
    amrex::MultiFab & acceleration,
    int level,
    amrex::Real timestep,
    amrex::Real a_old=1.0,
    amrex::Real a_half=1.0,
    int where_width=0
) 
```


Implements [*NyxParticleContainerBase::moveKickDrift*](classNyxParticleContainerBase.md#function-movekickdrift)


### <a href="#function-writeallatlevel" id="function-writeallatlevel">function writeAllAtLevel </a>


```cpp
void AGNParticleContainer::writeAllAtLevel (
    int lev
) 
```


Write out all particles at a level 


        

### <a href="#function-agnparticlecontainer" id="function-agnparticlecontainer">function ~AGNParticleContainer </a>


```cpp
inline virtual AGNParticleContainer::~AGNParticleContainer () 
```


## Protected Attributes Documentation


### <a href="#variable-real-comp-names" id="variable-real-comp-names">variable real\_comp\_names </a>


```cpp
amrex::Vector<std::string> AGNParticleContainer::real_comp_names;
```



### <a href="#variable-sub-cycle" id="variable-sub-cycle">variable sub\_cycle </a>


```cpp
bool AGNParticleContainer::sub_cycle;
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/AGNParticleContainer.H`