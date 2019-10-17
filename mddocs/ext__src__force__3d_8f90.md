
# File ext\_src\_force\_3d.f90


[**File List**](files.md) **>** [**Forcing**](dir_45682215f16eaf57f766b3c547de68bc.md) **>** [**ext\_src\_force\_3d.f90**](ext__src__force__3d_8f90.md)

[Go to the source code of this file.](ext__src__force__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**ext\_src\_force**](ext__src__force__3d_8f90.md#function-ext-src-force) (lo lo, hi hi, old\_state old\_state, os\_l1 os\_l1, os\_l2 os\_l2, os\_l3 os\_l3, os\_h1 os\_h1, os\_h2 os\_h2, os\_h3 os\_h3, new\_state new\_state, ns\_l1 ns\_l1, ns\_l2 ns\_l2, ns\_l3 ns\_l3, ns\_h1 ns\_h1, ns\_h2 ns\_h2, ns\_h3 ns\_h3, old\_diag old\_diag, od\_l1 od\_l1, od\_l2 od\_l2, od\_l3 od\_l3, od\_h1 od\_h1, od\_h2 od\_h2, od\_h3 od\_h3, new\_diag new\_diag, nd\_l1 nd\_l1, nd\_l2 nd\_l2, nd\_l3 nd\_l3, nd\_h1 nd\_h1, nd\_h2 nd\_h2, nd\_h3 nd\_h3, src src, src\_l1 src\_l1, src\_l2 src\_l2, src\_l3 src\_l3, src\_h1 src\_h1, src\_h2 src\_h2, src\_h3 src\_h3, problo problo, dx dx, time time, z z, dt dt) <br> |








## Public Functions Documentation


### <a href="#function-ext-src-force" id="function-ext-src-force">function ext\_src\_force </a>


```cpp
subroutine ext_src_force (
    lo lo,
    hi hi,
    old_state old_state,
    os_l1 os_l1,
    os_l2 os_l2,
    os_l3 os_l3,
    os_h1 os_h1,
    os_h2 os_h2,
    os_h3 os_h3,
    new_state new_state,
    ns_l1 ns_l1,
    ns_l2 ns_l2,
    ns_l3 ns_l3,
    ns_h1 ns_h1,
    ns_h2 ns_h2,
    ns_h3 ns_h3,
    old_diag old_diag,
    od_l1 od_l1,
    od_l2 od_l2,
    od_l3 od_l3,
    od_h1 od_h1,
    od_h2 od_h2,
    od_h3 od_h3,
    new_diag new_diag,
    nd_l1 nd_l1,
    nd_l2 nd_l2,
    nd_l3 nd_l3,
    nd_h1 nd_h1,
    nd_h2 nd_h2,
    nd_h3 nd_h3,
    src src,
    src_l1 src_l1,
    src_l2 src_l2,
    src_l3 src_l3,
    src_h1 src_h1,
    src_h2 src_h2,
    src_h3 src_h3,
    problo problo,
    dx dx,
    time time,
    z z,
    dt dt
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Forcing/ext_src_force_3d.f90`