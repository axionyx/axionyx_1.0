
# File Derive\_fdm\_3d.f90


[**File List**](files.md) **>** [**FDM**](dir_43b815edcf2a06ee60d8a45cc6c77fb8.md) **>** [**Derive\_fdm\_3d.f90**](Derive__fdm__3d_8f90.md)

[Go to the source code of this file.](Derive__fdm__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**ca\_axangmom\_x**](Derive__fdm__3d_8f90.md#function-ca-axangmom-x) (angmom\_x angmom\_x, angmom\_x\_l1 angmom\_x\_l1, angmom\_x\_l2 angmom\_x\_l2) <br> |
|  subroutine | [**ca\_axangmom\_y**](Derive__fdm__3d_8f90.md#function-ca-axangmom-y) (angmom\_y angmom\_y, angmom\_y\_l1 angmom\_y\_l1, angmom\_y\_l2 angmom\_y\_l2) <br> |
|  subroutine | [**ca\_axangmom\_z**](Derive__fdm__3d_8f90.md#function-ca-axangmom-z) (angmom\_z angmom\_z, angmom\_z\_l1 angmom\_z\_l1, angmom\_z\_l2 angmom\_z\_l2) <br> |
|  subroutine | [**ca\_axekin**](Derive__fdm__3d_8f90.md#function-ca-axekin) (ekin ekin, ekin\_l1 ekin\_l1, ekin\_l2 ekin\_l2, ekin\_l3 ekin\_l3, ekin\_h1 ekin\_h1, ekin\_h2 ekin\_h2) <br> |
|  subroutine | [**ca\_axekinrho**](Derive__fdm__3d_8f90.md#function-ca-axekinrho) (ekinrho ekinrho, ekinrho\_l1 ekinrho\_l1, ekinrho\_l2 ekinrho\_l2, ekinrho\_l3 ekinrho\_l3) <br> |
|  subroutine | [**ca\_axekinv**](Derive__fdm__3d_8f90.md#function-ca-axekinv) (ekinv ekinv, ekinv\_l1 ekinv\_l1, ekinv\_l2 ekinv\_l2, ekinv\_l3 ekinv\_l3, ekinv\_h1 ekinv\_h1) <br> |
|  subroutine | [**ca\_axepot**](Derive__fdm__3d_8f90.md#function-ca-axepot) (epot epot, epot\_l1 epot\_l1, epot\_l2 epot\_l2, epot\_l3 epot\_l3, epot\_h1 epot\_h1, epot\_h2 epot\_h2) <br> |
|  subroutine | [**ca\_axphase**](Derive__fdm__3d_8f90.md#function-ca-axphase) (phase phase, phase\_l1 phase\_l1, phase\_l2 phase\_l2, phase\_l3 phase\_l3, phase\_h1 phase\_h1) <br> |
|  subroutine | [**ca\_axvel**](Derive__fdm__3d_8f90.md#function-ca-axvel) (ekinv ekinv, ekinv\_l1 ekinv\_l1, ekinv\_l2 ekinv\_l2, ekinv\_l3 ekinv\_l3, ekinv\_h1 ekinv\_h1) <br> |
|  subroutine | [**ca\_dererrx**](Derive__fdm__3d_8f90.md#function-ca-dererrx) (err err, err\_l1 err\_l1, err\_l2 err\_l2, err\_l3 err\_l3, err\_h1 err\_h1, err\_h2 err\_h2) <br> |
|  subroutine | [**ca\_dererry**](Derive__fdm__3d_8f90.md#function-ca-dererry) (err err, err\_l1 err\_l1, err\_l2 err\_l2, err\_l3 err\_l3, err\_h1 err\_h1, err\_h2 err\_h2) <br> |
|  subroutine | [**ca\_dererrz**](Derive__fdm__3d_8f90.md#function-ca-dererrz) (err err, err\_l1 err\_l1, err\_l2 err\_l2, err\_l3 err\_l3, err\_h1 err\_h1, err\_h2 err\_h2) <br> |








## Public Functions Documentation


### <a href="#function-ca-axangmom-x" id="function-ca-axangmom-x">function ca\_axangmom\_x </a>


```cpp
subroutine ca_axangmom_x (
    angmom_x angmom_x,
    angmom_x_l1 angmom_x_l1,
    angmom_x_l2 angmom_x_l2
) 
```



### <a href="#function-ca-axangmom-y" id="function-ca-axangmom-y">function ca\_axangmom\_y </a>


```cpp
subroutine ca_axangmom_y (
    angmom_y angmom_y,
    angmom_y_l1 angmom_y_l1,
    angmom_y_l2 angmom_y_l2
) 
```



### <a href="#function-ca-axangmom-z" id="function-ca-axangmom-z">function ca\_axangmom\_z </a>


```cpp
subroutine ca_axangmom_z (
    angmom_z angmom_z,
    angmom_z_l1 angmom_z_l1,
    angmom_z_l2 angmom_z_l2
) 
```



### <a href="#function-ca-axekin" id="function-ca-axekin">function ca\_axekin </a>


```cpp
subroutine ca_axekin (
    ekin ekin,
    ekin_l1 ekin_l1,
    ekin_l2 ekin_l2,
    ekin_l3 ekin_l3,
    ekin_h1 ekin_h1,
    ekin_h2 ekin_h2
) 
```



### <a href="#function-ca-axekinrho" id="function-ca-axekinrho">function ca\_axekinrho </a>


```cpp
subroutine ca_axekinrho (
    ekinrho ekinrho,
    ekinrho_l1 ekinrho_l1,
    ekinrho_l2 ekinrho_l2,
    ekinrho_l3 ekinrho_l3
) 
```



### <a href="#function-ca-axekinv" id="function-ca-axekinv">function ca\_axekinv </a>


```cpp
subroutine ca_axekinv (
    ekinv ekinv,
    ekinv_l1 ekinv_l1,
    ekinv_l2 ekinv_l2,
    ekinv_l3 ekinv_l3,
    ekinv_h1 ekinv_h1
) 
```



### <a href="#function-ca-axepot" id="function-ca-axepot">function ca\_axepot </a>


```cpp
subroutine ca_axepot (
    epot epot,
    epot_l1 epot_l1,
    epot_l2 epot_l2,
    epot_l3 epot_l3,
    epot_h1 epot_h1,
    epot_h2 epot_h2
) 
```



### <a href="#function-ca-axphase" id="function-ca-axphase">function ca\_axphase </a>


```cpp
subroutine ca_axphase (
    phase phase,
    phase_l1 phase_l1,
    phase_l2 phase_l2,
    phase_l3 phase_l3,
    phase_h1 phase_h1
) 
```



### <a href="#function-ca-axvel" id="function-ca-axvel">function ca\_axvel </a>


```cpp
subroutine ca_axvel (
    ekinv ekinv,
    ekinv_l1 ekinv_l1,
    ekinv_l2 ekinv_l2,
    ekinv_l3 ekinv_l3,
    ekinv_h1 ekinv_h1
) 
```



### <a href="#function-ca-dererrx" id="function-ca-dererrx">function ca\_dererrx </a>


```cpp
subroutine ca_dererrx (
    err err,
    err_l1 err_l1,
    err_l2 err_l2,
    err_l3 err_l3,
    err_h1 err_h1,
    err_h2 err_h2
) 
```



### <a href="#function-ca-dererry" id="function-ca-dererry">function ca\_dererry </a>


```cpp
subroutine ca_dererry (
    err err,
    err_l1 err_l1,
    err_l2 err_l2,
    err_l3 err_l3,
    err_h1 err_h1,
    err_h2 err_h2
) 
```



### <a href="#function-ca-dererrz" id="function-ca-dererrz">function ca\_dererrz </a>


```cpp
subroutine ca_dererrz (
    err err,
    err_l1 err_l1,
    err_l2 err_l2,
    err_l3 err_l3,
    err_h1 err_h1,
    err_h2 err_h2
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/FDM/Derive_fdm_3d.f90`