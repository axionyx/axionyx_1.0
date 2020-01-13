
# File Nyx\_setup.cpp


[**File List**](files.md) **>** [**Initialization**](dir_71a4420ed1f8982e7234eb6a0b7e6d5d.md) **>** [**Nyx\_setup.cpp**](Nyx__setup_8cpp.md)

[Go to the source code of this file.](Nyx__setup_8cpp_source.md)



* `#include "AMReX_LevelBld.H"`
* `#include "Nyx.H"`
* `#include "Nyx_F.H"`
* `#include "Derive_F.H"`











## Public Types

| Type | Name |
| ---: | :--- |
| typedef StateDescriptor::BndryFunc | [**BndryFunc**](Nyx__setup_8cpp.md#typedef-bndryfunc)  <br> |



## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  int | [**norm\_vel\_bc**](Nyx__setup_8cpp.md#variable-norm-vel-bc)   = =
{
    INT\_DIR, EXT\_DIR, FOEXTRAP, REFLECT\_ODD, REFLECT\_ODD, REFLECT\_ODD
}<br> |
|  int | [**scalar\_bc**](Nyx__setup_8cpp.md#variable-scalar-bc)   = =
{
    INT\_DIR, EXT\_DIR, FOEXTRAP, REFLECT\_EVEN, REFLECT\_EVEN, REFLECT\_EVEN
}<br> |
|  int | [**tang\_vel\_bc**](Nyx__setup_8cpp.md#variable-tang-vel-bc)   = =
{
    INT\_DIR, EXT\_DIR, FOEXTRAP, REFLECT\_EVEN, REFLECT\_EVEN, REFLECT\_EVEN
}<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  Box | [**grow\_box\_by\_one**](Nyx__setup_8cpp.md#function-grow-box-by-one) (const Box & b) <br> |
|  Box | [**grow\_box\_by\_two**](Nyx__setup_8cpp.md#function-grow-box-by-two) (const Box & b) <br> |
|  void | [**set\_scalar\_bc**](Nyx__setup_8cpp.md#function-set-scalar-bc) (BCRec & bc, const BCRec & phys\_bc) <br> |
|  void | [**set\_x\_vel\_bc**](Nyx__setup_8cpp.md#function-set-x-vel-bc) (BCRec & bc, const BCRec & phys\_bc) <br> |
|  void | [**set\_y\_vel\_bc**](Nyx__setup_8cpp.md#function-set-y-vel-bc) (BCRec & bc, const BCRec & phys\_bc) <br> |
|  void | [**set\_z\_vel\_bc**](Nyx__setup_8cpp.md#function-set-z-vel-bc) (BCRec & bc, const BCRec & phys\_bc) <br> |
|  Box | [**the\_same\_box**](Nyx__setup_8cpp.md#function-the-same-box) (const Box & b) <br> |







## Public Types Documentation


### <a href="#typedef-bndryfunc" id="typedef-bndryfunc">typedef BndryFunc </a>


```cpp
typedef StateDescriptor::BndryFunc BndryFunc;
```


## Public Static Attributes Documentation


### <a href="#variable-norm-vel-bc" id="variable-norm-vel-bc">variable norm\_vel\_bc </a>


```cpp
int norm_vel_bc[];
```



### <a href="#variable-scalar-bc" id="variable-scalar-bc">variable scalar\_bc </a>


```cpp
int scalar_bc[];
```



### <a href="#variable-tang-vel-bc" id="variable-tang-vel-bc">variable tang\_vel\_bc </a>


```cpp
int tang_vel_bc[];
```


## Public Static Functions Documentation


### <a href="#function-grow-box-by-one" id="function-grow-box-by-one">function grow\_box\_by\_one </a>


```cpp
static Box grow_box_by_one (
    const Box & b
) 
```



### <a href="#function-grow-box-by-two" id="function-grow-box-by-two">function grow\_box\_by\_two </a>


```cpp
static Box grow_box_by_two (
    const Box & b
) 
```



### <a href="#function-set-scalar-bc" id="function-set-scalar-bc">function set\_scalar\_bc </a>


```cpp
static void set_scalar_bc (
    BCRec & bc,
    const BCRec & phys_bc
) 
```



### <a href="#function-set-x-vel-bc" id="function-set-x-vel-bc">function set\_x\_vel\_bc </a>


```cpp
static void set_x_vel_bc (
    BCRec & bc,
    const BCRec & phys_bc
) 
```



### <a href="#function-set-y-vel-bc" id="function-set-y-vel-bc">function set\_y\_vel\_bc </a>


```cpp
static void set_y_vel_bc (
    BCRec & bc,
    const BCRec & phys_bc
) 
```



### <a href="#function-set-z-vel-bc" id="function-set-z-vel-bc">function set\_z\_vel\_bc </a>


```cpp
static void set_z_vel_bc (
    BCRec & bc,
    const BCRec & phys_bc
) 
```



### <a href="#function-the-same-box" id="function-the-same-box">function the\_same\_box </a>


```cpp
static Box the_same_box (
    const Box & b
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Initialization/Nyx_setup.cpp`