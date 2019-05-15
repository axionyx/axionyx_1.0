
# File add\_grav\_source\_3d.f90


[**File List**](files.md) **>** [**HydroFortran**](dir_1fab266cd447ad3f3624320661f845f1.md) **>** [**add\_grav\_source\_3d.f90**](add__grav__source__3d_8f90.md)

[Go to the source code of this file.](add__grav__source__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**fort\_add\_grav\_source**](add__grav__source__3d_8f90.md#function-fort-add-grav-source) (lo lo, hi hi, uin uin, uin\_l1 uin\_l1, uin\_l2 uin\_l2, uin\_l3 uin\_l3, uin\_h1 uin\_h1, uin\_h2 uin\_h2, uin\_h3 uin\_h3, uout uout, uout\_l1 uout\_l1, uout\_l2 uout\_l2, uout\_l3 uout\_l3, uout\_h1 uout\_h1, uout\_h2 uout\_h2, uout\_h3 uout\_h3, grav grav, gv\_l1 gv\_l1, gv\_l2 gv\_l2, gv\_l3 gv\_l3, gv\_h1 gv\_h1, gv\_h2 gv\_h2, gv\_h3 gv\_h3, dt dt, a\_old a\_old, a\_new a\_new) <br> |








## Public Functions Documentation


### <a href="#function-fort-add-grav-source" id="function-fort-add-grav-source">function fort\_add\_grav\_source </a>


```cpp
subroutine fort_add_grav_source (
    lo lo,
    hi hi,
    uin uin,
    uin_l1 uin_l1,
    uin_l2 uin_l2,
    uin_l3 uin_l3,
    uin_h1 uin_h1,
    uin_h2 uin_h2,
    uin_h3 uin_h3,
    uout uout,
    uout_l1 uout_l1,
    uout_l2 uout_l2,
    uout_l3 uout_l3,
    uout_h1 uout_h1,
    uout_h2 uout_h2,
    uout_h3 uout_h3,
    grav grav,
    gv_l1 gv_l1,
    gv_l2 gv_l2,
    gv_l3 gv_l3,
    gv_h1 gv_h1,
    gv_h2 gv_h2,
    gv_h3 gv_h3,
    dt dt,
    a_old a_old,
    a_new a_new
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/rootfft/Source/HydroFortran/add_grav_source_3d.f90`