
# Class DarkMatterParticleContainer


[**Class List**](annotated.md) **>** [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md)








Inherits the following classes: [NyxParticleContainer](classNyxParticleContainer.md)








## Public Types

| Type | Name |
| ---: | :--- |
| typedef amrex::ParConstIter&lt; 1+BL\_SPACEDIM &gt; | [**MyConstParIter**](classDarkMatterParticleContainer.md#typedef-myconstpariter)  <br> |
| typedef amrex::ParIter&lt; 1+BL\_SPACEDIM &gt; | [**MyParIter**](classDarkMatterParticleContainer.md#typedef-mypariter)  <br> |

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
|  void | [**AssignDensityAndVels**](classDarkMatterParticleContainer.md#function-assigndensityandvels) (amrex::Vector&lt; std::unique\_ptr&lt; amrex::MultiFab &gt; &gt; & mf, int lev\_min=0) const<br> |
|   | [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md#function-darkmatterparticlecontainer) (amrex::Amr \* amr) <br> |
|  void | [**InitCosmo**](classDarkMatterParticleContainer.md#function-initcosmo-1-2) (amrex::MultiFab & mf, const amrex::Real vel\_fac, const amrex::Vector&lt; int &gt; n\_part, const amrex::Real particleMass) <br> |
|  void | [**InitCosmo**](classDarkMatterParticleContainer.md#function-initcosmo-2-2) (amrex::MultiFab & mf, const amrex::Real vel\_fac, const amrex::Vector&lt; int &gt; n\_part, const amrex::Real particleMass, const amrex::Real shift) <br> |
|  void | [**InitCosmo1ppc**](classDarkMatterParticleContainer.md#function-initcosmo1ppc) (amrex::MultiFab & mf, const amrex::Real vel\_fac, const amrex::Real particleMass) <br> |
|  void | [**InitCosmo1ppcMultiLevel**](classDarkMatterParticleContainer.md#function-initcosmo1ppcmultilevel) (amrex::MultiFab & mf, const amrex::Real disp\_fac, const amrex::Real vel\_fac, const amrex::Real particleMass, int disp\_idx, int vel\_idx, amrex::BoxArray & baWhereNot, int lev, int nlevs) <br> |
|  void | [**InitFromBinaryMortonFile**](classDarkMatterParticleContainer.md#function-initfrombinarymortonfile) (const std::string & particle\_directory, int nextra, int skip\_factor) <br> |
| virtual void | [**moveKick**](classDarkMatterParticleContainer.md#function-movekick) (amrex::MultiFab & acceleration, int level, amrex::Real timestep, amrex::Real a\_new=1.0, amrex::Real a\_half=1.0) <br> |
| virtual void | [**moveKickDrift**](classDarkMatterParticleContainer.md#function-movekickdrift) (amrex::MultiFab & acceleration, int level, amrex::Real timestep, amrex::Real a\_old=1.0, amrex::Real a\_half=1.0, int where\_width=0) <br> |
| virtual  | [**~DarkMatterParticleContainer**](classDarkMatterParticleContainer.md#function-darkmatterparticlecontainer) () <br> |

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











## Protected Attributes inherited from [NyxParticleContainer](classNyxParticleContainer.md)

| Type | Name |
| ---: | :--- |
|  amrex::Vector&lt; std::string &gt; | [**real\_comp\_names**](classNyxParticleContainer.md#variable-real-comp-names)  <br> |
|  bool | [**sub\_cycle**](classNyxParticleContainer.md#variable-sub-cycle)  <br> |











## Public Types Documentation


### <a href="#typedef-myconstpariter" id="typedef-myconstpariter">typedef MyConstParIter </a>


```cpp
using DarkMatterParticleContainer::MyConstParIter =  amrex::ParConstIter<1+BL_SPACEDIM>;
```



### <a href="#typedef-mypariter" id="typedef-mypariter">typedef MyParIter </a>


```cpp
using DarkMatterParticleContainer::MyParIter =  amrex::ParIter<1+BL_SPACEDIM>;
```


## Public Functions Documentation


### <a href="#function-assigndensityandvels" id="function-assigndensityandvels">function AssignDensityAndVels </a>


```cpp
void DarkMatterParticleContainer::AssignDensityAndVels (
    amrex::Vector< std::unique_ptr< amrex::MultiFab > > & mf,
    int lev_min=0
) const
```



### <a href="#function-darkmatterparticlecontainer" id="function-darkmatterparticlecontainer">function DarkMatterParticleContainer </a>


```cpp
inline DarkMatterParticleContainer::DarkMatterParticleContainer (
    amrex::Amr * amr
) 
```



### <a href="#function-initcosmo-1-2" id="function-initcosmo-1-2">function InitCosmo [1/2]</a>


```cpp
void DarkMatterParticleContainer::InitCosmo (
    amrex::MultiFab & mf,
    const amrex::Real vel_fac,
    const amrex::Vector< int > n_part,
    const amrex::Real particleMass
) 
```



### <a href="#function-initcosmo-2-2" id="function-initcosmo-2-2">function InitCosmo [2/2]</a>


```cpp
void DarkMatterParticleContainer::InitCosmo (
    amrex::MultiFab & mf,
    const amrex::Real vel_fac,
    const amrex::Vector< int > n_part,
    const amrex::Real particleMass,
    const amrex::Real shift
) 
```



### <a href="#function-initcosmo1ppc" id="function-initcosmo1ppc">function InitCosmo1ppc </a>


```cpp
void DarkMatterParticleContainer::InitCosmo1ppc (
    amrex::MultiFab & mf,
    const amrex::Real vel_fac,
    const amrex::Real particleMass
) 
```



### <a href="#function-initcosmo1ppcmultilevel" id="function-initcosmo1ppcmultilevel">function InitCosmo1ppcMultiLevel </a>


```cpp
void DarkMatterParticleContainer::InitCosmo1ppcMultiLevel (
    amrex::MultiFab & mf,
    const amrex::Real disp_fac,
    const amrex::Real vel_fac,
    const amrex::Real particleMass,
    int disp_idx,
    int vel_idx,
    amrex::BoxArray & baWhereNot,
    int lev,
    int nlevs
) 
```



### <a href="#function-initfrombinarymortonfile" id="function-initfrombinarymortonfile">function InitFromBinaryMortonFile </a>


```cpp
void DarkMatterParticleContainer::InitFromBinaryMortonFile (
    const std::string & particle_directory,
    int nextra,
    int skip_factor
) 
```



### <a href="#function-movekick" id="function-movekick">function moveKick </a>


```cpp
virtual void DarkMatterParticleContainer::moveKick (
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
virtual void DarkMatterParticleContainer::moveKickDrift (
    amrex::MultiFab & acceleration,
    int level,
    amrex::Real timestep,
    amrex::Real a_old=1.0,
    amrex::Real a_half=1.0,
    int where_width=0
) 
```


Implements [*NyxParticleContainerBase::moveKickDrift*](classNyxParticleContainerBase.md#function-movekickdrift)


### <a href="#function-darkmatterparticlecontainer" id="function-darkmatterparticlecontainer">function ~DarkMatterParticleContainer </a>


```cpp
inline virtual DarkMatterParticleContainer::~DarkMatterParticleContainer () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/DarkMatterParticleContainer.H`