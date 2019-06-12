
# Namespace anonymous\_namespace{NyxParticles.cpp}


[**Class List**](annotated.md) **>** [**anonymous\_namespace{NyxParticles.cpp}**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md)


















## Public Attributes

| Type | Name |
| ---: | :--- |
|  Vector&lt; [**NyxParticleContainerBase**](classNyxParticleContainerBase.md) \* &gt; | [**ActiveParticles**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-activeparticles)  <br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**DMPC**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-dmpc)   = = 0<br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**GhostPC**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-ghostpc)   = = 0<br> |
|  Vector&lt; [**NyxParticleContainerBase**](classNyxParticleContainerBase.md) \* &gt; | [**GhostParticles**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-ghostparticles)  <br> |
|  [**StellarParticleContainer**](classNyxParticleContainer.md) \* | [**GhostSPC**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-ghostspc)   = = 0<br> |
|  [**StellarParticleContainer**](classNyxParticleContainer.md) \* | [**SPC**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-spc)   = = 0<br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**SPHPC**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-sphpc)   = = 0<br> |
|  [**DarkMatterParticleContainer**](classDarkMatterParticleContainer.md) \* | [**VirtPC**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-virtpc)   = = 0<br> |
|  [**StellarParticleContainer**](classNyxParticleContainer.md) \* | [**VirtSPC**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-virtspc)   = = 0<br> |
|  Vector&lt; [**NyxParticleContainerBase**](classNyxParticleContainerBase.md) \* &gt; | [**VirtualParticles**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-virtualparticles)  <br> |
|  std::string | [**ascii\_particle\_file**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-ascii-particle-file)  <br> |
|  std::string | [**binary\_particle\_file**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-binary-particle-file)  <br> |
|  std::string | [**sph\_particle\_file**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-sph-particle-file)  <br> |
|  bool | [**virtual\_particles\_set**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#variable-virtual-particles-set)   = = false<br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**RemoveParticlesOnExit**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#function-removeparticlesonexit) () <br> |
|  const std::string | [**agn\_chk\_particle\_file**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#function-agn-chk-particle-file) ("AGN") <br> |
|  const std::string | [**dm\_chk\_particle\_file**](namespaceanonymous__namespace_02NyxParticles_8cpp_03.md#function-dm-chk-particle-file) ("DM") <br> |








## Public Attributes Documentation


### <a href="#variable-activeparticles" id="variable-activeparticles">variable ActiveParticles </a>


```cpp
Vector<NyxParticleContainerBase*> anonymous_namespace{NyxParticles.cpp}::ActiveParticles;
```



### <a href="#variable-dmpc" id="variable-dmpc">variable DMPC </a>


```cpp
DarkMatterParticleContainer* anonymous_namespace{NyxParticles.cpp}::DMPC;
```



### <a href="#variable-ghostpc" id="variable-ghostpc">variable GhostPC </a>


```cpp
DarkMatterParticleContainer* anonymous_namespace{NyxParticles.cpp}::GhostPC;
```



### <a href="#variable-ghostparticles" id="variable-ghostparticles">variable GhostParticles </a>


```cpp
Vector<NyxParticleContainerBase*> anonymous_namespace{NyxParticles.cpp}::GhostParticles;
```



### <a href="#variable-ghostspc" id="variable-ghostspc">variable GhostSPC </a>


```cpp
StellarParticleContainer* anonymous_namespace{NyxParticles.cpp}::GhostSPC;
```



### <a href="#variable-spc" id="variable-spc">variable SPC </a>


```cpp
StellarParticleContainer* anonymous_namespace{NyxParticles.cpp}::SPC;
```



### <a href="#variable-sphpc" id="variable-sphpc">variable SPHPC </a>


```cpp
DarkMatterParticleContainer* anonymous_namespace{NyxParticles.cpp}::SPHPC;
```



### <a href="#variable-virtpc" id="variable-virtpc">variable VirtPC </a>


```cpp
DarkMatterParticleContainer* anonymous_namespace{NyxParticles.cpp}::VirtPC;
```



### <a href="#variable-virtspc" id="variable-virtspc">variable VirtSPC </a>


```cpp
StellarParticleContainer* anonymous_namespace{NyxParticles.cpp}::VirtSPC;
```



### <a href="#variable-virtualparticles" id="variable-virtualparticles">variable VirtualParticles </a>


```cpp
Vector<NyxParticleContainerBase*> anonymous_namespace{NyxParticles.cpp}::VirtualParticles;
```



### <a href="#variable-ascii-particle-file" id="variable-ascii-particle-file">variable ascii\_particle\_file </a>


```cpp
std::string anonymous_namespace{NyxParticles.cpp}::ascii_particle_file;
```



### <a href="#variable-binary-particle-file" id="variable-binary-particle-file">variable binary\_particle\_file </a>


```cpp
std::string anonymous_namespace{NyxParticles.cpp}::binary_particle_file;
```



### <a href="#variable-sph-particle-file" id="variable-sph-particle-file">variable sph\_particle\_file </a>


```cpp
std::string anonymous_namespace{NyxParticles.cpp}::sph_particle_file;
```



### <a href="#variable-virtual-particles-set" id="variable-virtual-particles-set">variable virtual\_particles\_set </a>


```cpp
bool anonymous_namespace{NyxParticles.cpp}::virtual_particles_set;
```


## Public Functions Documentation


### <a href="#function-removeparticlesonexit" id="function-removeparticlesonexit">function RemoveParticlesOnExit </a>


```cpp
void anonymous_namespace{NyxParticles.cpp}::RemoveParticlesOnExit () 
```



### <a href="#function-agn-chk-particle-file" id="function-agn-chk-particle-file">function agn\_chk\_particle\_file </a>


```cpp
const std::string anonymous_namespace{NyxParticles.cpp}::agn_chk_particle_file (
    "AGN"
) 
```



### <a href="#function-dm-chk-particle-file" id="function-dm-chk-particle-file">function dm\_chk\_particle\_file </a>


```cpp
const std::string anonymous_namespace{NyxParticles.cpp}::dm_chk_particle_file (
    "DM"
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/NyxParticles.cpp`