
# File f\_rhs.f90


[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**f\_rhs.f90**](f__rhs_8f90.md)

[Go to the source code of this file.](f__rhs_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**f\_rhs**](f__rhs_8f90.md#function-f-rhs) (num\_eq num\_eq, time time, e\_in e\_in, energy energy, rpar rpar, ipar ipar) <br> |
|  subroutine | [**f\_rhs\_rpar**](f__rhs_8f90.md#function-f-rhs-rpar) (num\_eq num\_eq, time time, e\_in e\_in, energy energy, rpar rpar, ipar ipar) <br> |
|  subroutine | [**f\_rhs\_split**](f__rhs_8f90.md#function-f-rhs-split) (num\_eq num\_eq, time time, y\_in y\_in, yp\_out yp\_out, rpar rpar, ipar ipar) <br> |
|  subroutine | [**f\_rhs\_vec**](f__rhs_8f90.md#function-f-rhs-vec) (time time, e\_in e\_in, energy energy) <br> |
|  subroutine | [**jac**](f__rhs_8f90.md#function-jac) (neq neq, t t, y y, ml ml, mu mu, pd pd, nrpd nrpd, rpar rpar, ipar ipar) <br> |








## Public Functions Documentation


### <a href="#function-f-rhs" id="function-f-rhs">function f\_rhs </a>


```cpp
subroutine f_rhs (
    num_eq num_eq,
    time time,
    e_in e_in,
    energy energy,
    rpar rpar,
    ipar ipar
) 
```



### <a href="#function-f-rhs-rpar" id="function-f-rhs-rpar">function f\_rhs\_rpar </a>


```cpp
subroutine f_rhs_rpar (
    num_eq num_eq,
    time time,
    e_in e_in,
    energy energy,
    rpar rpar,
    ipar ipar
) 
```



### <a href="#function-f-rhs-split" id="function-f-rhs-split">function f\_rhs\_split </a>


```cpp
subroutine f_rhs_split (
    num_eq num_eq,
    time time,
    y_in y_in,
    yp_out yp_out,
    rpar rpar,
    ipar ipar
) 
```



### <a href="#function-f-rhs-vec" id="function-f-rhs-vec">function f\_rhs\_vec </a>


```cpp
subroutine f_rhs_vec (
    time time,
    e_in e_in,
    energy energy
) 
```



### <a href="#function-jac" id="function-jac">function jac </a>


```cpp
subroutine jac (
    neq neq,
    t t,
    y y,
    ml ml,
    mu mu,
    pd pd,
    nrpd nrpd,
    rpar rpar,
    ipar ipar
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/HeatCool/f_rhs.f90`