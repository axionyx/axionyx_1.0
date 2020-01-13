
# File Tagging\_3d.f90


[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Tagging**](dir_c14a965952b26c2f69053cc66c8fb69f.md) **>** [**Tagging\_3d.f90**](Source_2Tagging_2Tagging__3d_8f90.md)

[Go to the source code of this file.](Source_2Tagging_2Tagging__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**tag\_denerror**](Source_2Tagging_2Tagging__3d_8f90.md#function-tag-denerror) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2, tagh3 tagh3, set set, clear clear, den den, denl1 denl1, denl2 denl2, denl3 denl3, denh1 denh1, denh2 denh2, denh3 denh3, lo lo, hi hi, nd nd, domlo domlo, domhi domhi, delta delta, xlo xlo, problo problo, time time, level level) <br> |
|  subroutine | [**tag\_laplac\_error**](Source_2Tagging_2Tagging__3d_8f90.md#function-tag-laplac-error) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2) <br> |
|  subroutine | [**tag\_overdensity**](Source_2Tagging_2Tagging__3d_8f90.md#function-tag-overdensity) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2) <br> |
|  subroutine | [**tag\_part\_cnt\_err**](Source_2Tagging_2Tagging__3d_8f90.md#function-tag-part-cnt-err) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2) <br> |
|  subroutine | [**tag\_presserror**](Source_2Tagging_2Tagging__3d_8f90.md#function-tag-presserror) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2, tagh3 tagh3) <br> |
|  subroutine | [**tag\_region**](Source_2Tagging_2Tagging__3d_8f90.md#function-tag-region) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2, tagh3 tagh3, set set, lo lo, hi hi, domlo domlo, domhi domhi, delta delta, xlo xlo, problo problo, level level) <br> |
|  subroutine | [**tag\_temperror**](Source_2Tagging_2Tagging__3d_8f90.md#function-tag-temperror) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2, tagh3 tagh3) <br> |
|  subroutine | [**tag\_velerror**](Source_2Tagging_2Tagging__3d_8f90.md#function-tag-velerror) (tag tag, tagl1 tagl1, tagl2 tagl2, tagl3 tagl3, tagh1 tagh1, tagh2 tagh2, tagh3 tagh3, set set, clear clear, vel vel, vell1 vell1, vell2 vell2, vell3 vell3, velh1 velh1, velh2 velh2, velh3 velh3, lo lo, hi hi, nv nv, domlo domlo, domhi domhi, delta delta, xlo xlo, problo problo, time time, level level) <br> |








## Public Functions Documentation


### <a href="#function-tag-denerror" id="function-tag-denerror">function tag\_denerror </a>


```cpp
subroutine tag_denerror (
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
    nd nd,
    domlo domlo,
    domhi domhi,
    delta delta,
    xlo xlo,
    problo problo,
    time time,
    level level
) 
```



### <a href="#function-tag-laplac-error" id="function-tag-laplac-error">function tag\_laplac\_error </a>


```cpp
subroutine tag_laplac_error (
    tag tag,
    tagl1 tagl1,
    tagl2 tagl2,
    tagl3 tagl3,
    tagh1 tagh1,
    tagh2 tagh2
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



### <a href="#function-tag-part-cnt-err" id="function-tag-part-cnt-err">function tag\_part\_cnt\_err </a>


```cpp
subroutine tag_part_cnt_err (
    tag tag,
    tagl1 tagl1,
    tagl2 tagl2,
    tagl3 tagl3,
    tagh1 tagh1,
    tagh2 tagh2
) 
```



### <a href="#function-tag-presserror" id="function-tag-presserror">function tag\_presserror </a>


```cpp
subroutine tag_presserror (
    tag tag,
    tagl1 tagl1,
    tagl2 tagl2,
    tagl3 tagl3,
    tagh1 tagh1,
    tagh2 tagh2,
    tagh3 tagh3
) 
```



### <a href="#function-tag-region" id="function-tag-region">function tag\_region </a>


```cpp
subroutine tag_region (
    tag tag,
    tagl1 tagl1,
    tagl2 tagl2,
    tagl3 tagl3,
    tagh1 tagh1,
    tagh2 tagh2,
    tagh3 tagh3,
    set set,
    lo lo,
    hi hi,
    domlo domlo,
    domhi domhi,
    delta delta,
    xlo xlo,
    problo problo,
    level level
) 
```



### <a href="#function-tag-temperror" id="function-tag-temperror">function tag\_temperror </a>


```cpp
subroutine tag_temperror (
    tag tag,
    tagl1 tagl1,
    tagl2 tagl2,
    tagl3 tagl3,
    tagh1 tagh1,
    tagh2 tagh2,
    tagh3 tagh3
) 
```



### <a href="#function-tag-velerror" id="function-tag-velerror">function tag\_velerror </a>


```cpp
subroutine tag_velerror (
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
    nv nv,
    domlo domlo,
    domhi domhi,
    delta delta,
    xlo xlo,
    problo problo,
    time time,
    level level
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Tagging/Tagging_3d.f90`