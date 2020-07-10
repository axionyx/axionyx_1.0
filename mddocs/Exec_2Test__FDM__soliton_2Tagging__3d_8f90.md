
# File Tagging\_3d.f90


[**File List**](files.md) **>** [**Exec**](dir_43a12cefb7942b6f49b5b628aafd3192.md) **>** [**Test\_FDM\_soliton**](dir_5cdd3d89825e7a38a6ef6840955c84ea.md) **>** [**Tagging\_3d.f90**](Exec_2Test__FDM__soliton_2Tagging__3d_8f90.md)

[Go to the source code of this file.](Exec_2Test__FDM__soliton_2Tagging__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**tag\_axvel**](Exec_2Test__FDM__soliton_2Tagging__3d_8f90.md#function-tag-axvel) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2, tagh3 tagh3, set set, clear clear, vel vel, vell1 vell1, vell2 vell2, vell3 vell3, velh1 velh1, velh2 velh2, velh3 velh3, lo lo, hi hi, nc nc, domlo domlo, domhi domhi, delta delta, xlo xlo, problo problo) <br> |
|  subroutine | [**tag\_center**](Exec_2Test__FDM__soliton_2Tagging__3d_8f90.md#function-tag-center) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2, tagh3 tagh3, set set, clear clear, den den, denl1 denl1, denl2 denl2, denl3 denl3, denh1 denh1, denh2 denh2, denh3 denh3, lo lo, hi hi, nc nc, domlo domlo, domhi domhi, delta delta, level level, avg\_den avg\_den) <br> |
|  subroutine | [**tag\_overdensity**](Exec_2Test__FDM__soliton_2Tagging__3d_8f90.md#function-tag-overdensity) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2) <br> |








## Public Functions Documentation


### <a href="#function-tag-axvel" id="function-tag-axvel">function tag\_axvel </a>


```cpp
subroutine tag_axvel (
    tag tag,
    tagl1 tagl1,
    tagl2 tagl2,
    tagl3 tagl3,
    tagh1 tagh1,
    tagh2 tagh2,
    tagh3 tagh3,
    set set,
    clear clear,
    vel vel,
    vell1 vell1,
    vell2 vell2,
    vell3 vell3,
    velh1 velh1,
    velh2 velh2,
    velh3 velh3,
    lo lo,
    hi hi,
    nc nc,
    domlo domlo,
    domhi domhi,
    delta delta,
    xlo xlo,
    problo problo
) 
```



### <a href="#function-tag-center" id="function-tag-center">function tag\_center </a>


```cpp
subroutine tag_center (
    tag tag,
    tagl1 tagl1,
    tagl2 tagl2,
    tagl3 tagl3,
    tagh1 tagh1,
    tagh2 tagh2,
    tagh3 tagh3,
    set set,
    clear clear,
    den den,
    denl1 denl1,
    denl2 denl2,
    denl3 denl3,
    denh1 denh1,
    denh2 denh2,
    denh3 denh3,
    lo lo,
    hi hi,
    nc nc,
    domlo domlo,
    domhi domhi,
    delta delta,
    level level,
    avg_den avg_den
) 
```



### <a href="#function-tag-overdensity" id="function-tag-overdensity">function tag\_overdensity </a>


```cpp
subroutine tag_overdensity (
    tag tag,
    tagl1 tagl1,
    tagl2 tagl2,
    tagl3 tagl3,
    tagh1 tagh1,
    tagh2 tagh2
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Exec/Test_FDM_soliton/Tagging_3d.f90`