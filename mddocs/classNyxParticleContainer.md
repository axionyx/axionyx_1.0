
# Class NyxParticleContainer

**template &lt;int NSR, int NSI, int NAR, int NAI&gt;**


[**Class List**](annotated.md) **>** [**NyxParticleContainer**](classNyxParticleContainer.md)








Inherits the following classes: amrex::NeighborParticleContainer< NSR, NSI >,  [NyxParticleContainerBase](classNyxParticleContainerBase.md)








## Public Types

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
|  amrex::Vector&lt; std::string &gt; | [**real\_comp\_names**](classNyxParticleContainer.md#variable-real-comp-names)  <br> |
|  bool | [**sub\_cycle**](classNyxParticleContainer.md#variable-sub-cycle)  <br> |








## Public Types Documentation


### <a href="#typedef-aos" id="typedef-aos">typedef AoS </a>


```cpp
using NyxParticleContainer< NSR, NSI, NAR, NAI >::AoS =  typename amrex::ParticleContainer<NSR,NSI,NAR,NAI>::AoS;
```



### <a href="#typedef-myconstpariter" id="typedef-myconstpariter">typedef MyConstParIter </a>


```cpp
using NyxParticleContainer< NSR, NSI, NAR, NAI >::MyConstParIter =  amrex::ParConstIter<NSR,NSI,NAR,NAI>;
```



### <a href="#typedef-mypariter" id="typedef-mypariter">typedef MyParIter </a>


```cpp
using NyxParticleContainer< NSR, NSI, NAR, NAI >::MyParIter =  amrex::ParIter<NSR,NSI,NAR,NAI>;
```



### <a href="#typedef-particlelevel" id="typedef-particlelevel">typedef ParticleLevel </a>


```cpp
using NyxParticleContainer< NSR, NSI, NAR, NAI >::ParticleLevel =  typename amrex::ParticleContainer<NSR,NSI,NAR,NAI>::ParticleLevel;
```



### <a href="#typedef-particletiletype" id="typedef-particletiletype">typedef ParticleTileType </a>


```cpp
using NyxParticleContainer< NSR, NSI, NAR, NAI >::ParticleTileType =  amrex::ParticleTile<NSR,NSI,NAR,NAI>;
```



### <a href="#typedef-particletype" id="typedef-particletype">typedef ParticleType </a>


```cpp
typedef amrex::Particle<NSR,NSI> NyxParticleContainer< NSR, NSI, NAR, NAI >::ParticleType;
```


## Public Functions Documentation


### <a href="#function-assigndensity" id="function-assigndensity">function AssignDensity </a>


```cpp
inline virtual void NyxParticleContainer::AssignDensity (
    amrex::Vector< std::unique_ptr< amrex::MultiFab > > & mf,
    int lev_min=0,
    int ncomp=1,
    int finest_level=-1,
    int ngrow=1
) override const
```


Implements [*NyxParticleContainerBase::AssignDensity*](classNyxParticleContainerBase.md#function-assigndensity)


### <a href="#function-assigndensitysinglelevel" id="function-assigndensitysinglelevel">function AssignDensitySingleLevel </a>


```cpp
inline virtual void NyxParticleContainer::AssignDensitySingleLevel (
    amrex::MultiFab & mf,
    int level,
    int ncomp=1,
    int particle_lvl_offset=0
) override const
```


Implements [*NyxParticleContainerBase::AssignDensitySingleLevel*](classNyxParticleContainerBase.md#function-assigndensitysinglelevel)


### <a href="#function-getparticlevelocities" id="function-getparticlevelocities">function GetParticleVelocities </a>


```cpp
void NyxParticleContainer::GetParticleVelocities (
    amrex::Vector< amrex::Real > & part_vels
) 
```



### <a href="#function-multiplyparticlemass" id="function-multiplyparticlemass">function MultiplyParticleMass </a>


```cpp
void NyxParticleContainer::MultiplyParticleMass (
    int lev,
    amrex::Real mult
) 
```



### <a href="#function-nyxcheckpoint" id="function-nyxcheckpoint">function NyxCheckpoint </a>


```cpp
virtual void NyxParticleContainer::NyxCheckpoint (
    const std::string & dir,
    const std::string & name
) const
```



### <a href="#function-nyxparticlecontainer" id="function-nyxparticlecontainer">function NyxParticleContainer </a>


```cpp
inline NyxParticleContainer::NyxParticleContainer (
    amrex::Amr * amr,
    int nghost=0
) 
```



### <a href="#function-redistribute" id="function-redistribute">function Redistribute </a>


```cpp
inline virtual void NyxParticleContainer::Redistribute (
    int lev_min=0,
    int lev_max=-1,
    int nGrow=0
) override
```


Implements [*NyxParticleContainerBase::Redistribute*](classNyxParticleContainerBase.md#function-redistribute)


### <a href="#function-removeparticlesatlevel" id="function-removeparticlesatlevel">function RemoveParticlesAtLevel </a>


```cpp
inline virtual void NyxParticleContainer::RemoveParticlesAtLevel (
    int level
) override
```


Implements [*NyxParticleContainerBase::RemoveParticlesAtLevel*](classNyxParticleContainerBase.md#function-removeparticlesatlevel)


### <a href="#function-setparticlevelocities" id="function-setparticlevelocities">function SetParticleVelocities </a>


```cpp
void NyxParticleContainer::SetParticleVelocities (
    amrex::Vector< amrex::Real > & part_data
) 
```



### <a href="#function-writenyxplotfile" id="function-writenyxplotfile">function WriteNyxPlotFile </a>


```cpp
virtual void NyxParticleContainer::WriteNyxPlotFile (
    const std::string & dir,
    const std::string & name
) const
```



### <a href="#function-esttimestep-1-2" id="function-esttimestep-1-2">function estTimestep [1/2]</a>


```cpp
amrex::Real NyxParticleContainer::estTimestep (
    amrex::MultiFab & acceleration,
    int level,
    amrex::Real cfl
) const
```



### <a href="#function-esttimestep-2-2" id="function-esttimestep-2-2">function estTimestep [2/2]</a>


```cpp
amrex::Real NyxParticleContainer::estTimestep (
    amrex::MultiFab & acceleration,
    amrex::Real a,
    int level,
    amrex::Real cfl
) const
```



### <a href="#function-finestlevel" id="function-finestlevel">function finestLevel </a>


```cpp
inline virtual int NyxParticleContainer::finestLevel () override const
```


Implements [*NyxParticleContainerBase::finestLevel*](classNyxParticleContainerBase.md#function-finestlevel)


### <a href="#function-sumparticlemass" id="function-sumparticlemass">function sumParticleMass </a>


```cpp
inline virtual amrex::Real NyxParticleContainer::sumParticleMass (
    int level
) override const
```


Implements [*NyxParticleContainerBase::sumParticleMass*](classNyxParticleContainerBase.md#function-sumparticlemass)


### <a href="#function-sumparticlemomentum" id="function-sumparticlemomentum">function sumParticleMomentum </a>


```cpp
void NyxParticleContainer::sumParticleMomentum (
    int lev,
    amrex::Real * mom
) const
```



### <a href="#function-nyxparticlecontainer" id="function-nyxparticlecontainer">function ~NyxParticleContainer </a>


```cpp
inline virtual NyxParticleContainer::~NyxParticleContainer () 
```


## Protected Attributes Documentation


### <a href="#variable-real-comp-names" id="variable-real-comp-names">variable real\_comp\_names </a>


```cpp
amrex::Vector<std::string> NyxParticleContainer< NSR, NSI, NAR, NAI >::real_comp_names;
```



### <a href="#variable-sub-cycle" id="variable-sub-cycle">variable sub\_cycle </a>


```cpp
bool NyxParticleContainer< NSR, NSI, NAR, NAI >::sub_cycle;
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/NyxParticleContainer.H`