
# Namespace anonymous\_namespace{DarkMatterParticleContainer.cpp}


[**Class List**](annotated.md) **>** [**anonymous\_namespace{DarkMatterParticleContainer.cpp}**](namespaceanonymous__namespace_02DarkMatterParticleContainer_8cpp_03.md)



[More...](#detailed-description)











## Classes

| Type | Name |
| ---: | :--- |
| struct | [**BoxMortonKey**](structanonymous__namespace_02DarkMatterParticleContainer_8cpp_03_1_1BoxMortonKey.md) <br> |
| struct | [**ParticleMortonFileHeader**](structanonymous__namespace_02DarkMatterParticleContainer_8cpp_03_1_1ParticleMortonFileHeader.md) <br> |
| struct | [**by\_morton\_id**](structanonymous__namespace_02DarkMatterParticleContainer_8cpp_03_1_1by__morton__id.md) <br> |





## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**ReadHeader**](namespaceanonymous__namespace_02DarkMatterParticleContainer_8cpp_03.md#function-readheader) (const std::string & dir, const std::string & file, [**ParticleMortonFileHeader**](structanonymous__namespace_02DarkMatterParticleContainer_8cpp_03_1_1ParticleMortonFileHeader.md) & hdr) <br> |
|  std::string | [**get\_file\_name**](namespaceanonymous__namespace_02DarkMatterParticleContainer_8cpp_03.md#function-get-file-name) (const std::string & base, int file\_num) <br> |
|  uint64\_t | [**get\_morton\_index**](namespaceanonymous__namespace_02DarkMatterParticleContainer_8cpp_03.md#function-get-morton-index) (unsigned int x, unsigned int y, unsigned int z) <br> |
|  uint64\_t | [**split**](namespaceanonymous__namespace_02DarkMatterParticleContainer_8cpp_03.md#function-split) (unsigned int a) <br> |








## Public Functions Documentation


### <a href="#function-readheader" id="function-readheader">function ReadHeader </a>


```cpp
void anonymous_namespace{DarkMatterParticleContainer.cpp}::ReadHeader (
    const std::string & dir,
    const std::string & file,
    ParticleMortonFileHeader & hdr
) 
```



### <a href="#function-get-file-name" id="function-get-file-name">function get\_file\_name </a>


```cpp
std::string anonymous_namespace{DarkMatterParticleContainer.cpp}::get_file_name (
    const std::string & base,
    int file_num
) 
```



### <a href="#function-get-morton-index" id="function-get-morton-index">function get\_morton\_index </a>


```cpp
inline uint64_t anonymous_namespace{DarkMatterParticleContainer.cpp}::get_morton_index (
    unsigned int x,
    unsigned int y,
    unsigned int z
) 
```



### <a href="#function-split" id="function-split">function split </a>


```cpp
inline uint64_t anonymous_namespace{DarkMatterParticleContainer.cpp}::split (
    unsigned int a
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/DarkMatterParticleContainer.cpp`