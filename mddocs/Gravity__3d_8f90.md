
# File Gravity\_3d.f90


[**File List**](files.md) **>** [**Gravity**](dir_fdbf5007869eac89a42b1cd44aeda050.md) **>** [**Gravity\_3d.f90**](Gravity__3d_8f90.md)

[Go to the source code of this file.](Gravity__3d_8f90_source.md)


















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**fort\_avgdown\_phi**](Gravity__3d_8f90.md#function-fort-avgdown-phi) (crse crse, c\_l1 c\_l1, c\_l2 c\_l2, c\_l3 c\_l3, c\_h1 c\_h1, c\_h2 c\_h2, c\_h3 c\_h3, fine fine, f\_l1 f\_l1, f\_l2 f\_l2, f\_l3 f\_l3, f\_h1 f\_h1, f\_h2 f\_h2, f\_h3 f\_h3, lo lo, hi hi, lrat lrat) <br> |
|  subroutine | [**fort\_edge\_interp**](Gravity__3d_8f90.md#function-fort-edge-interp) (flo flo, fhi fhi, nc nc, ratio ratio, dir dir, fine fine, fine\_l0 fine\_l0, fine\_l1 fine\_l1, fine\_l2 fine\_l2, fine\_h0 fine\_h0, fine\_h1 fine\_h1, fine\_h2 fine\_h2) <br> |
|  subroutine | [**fort\_pc\_edge\_interp**](Gravity__3d_8f90.md#function-fort-pc-edge-interp) (lo lo, hi hi, nc nc, ratio ratio, dir dir, crse crse, crse\_l0 crse\_l0, crse\_l1 crse\_l1, crse\_l2 crse\_l2, crse\_h0 crse\_h0, crse\_h1 crse\_h1, crse\_h2 crse\_h2, fine fine, fine\_l0 fine\_l0, fine\_l1 fine\_l1, fine\_l2 fine\_l2, fine\_h0 fine\_h0, fine\_h1 fine\_h1, fine\_h2 fine\_h2) <br> |








## Public Functions Documentation


### <a href="#function-fort-avgdown-phi" id="function-fort-avgdown-phi">function fort\_avgdown\_phi </a>


```cpp
subroutine fort_avgdown_phi (
    crse crse,
    c_l1 c_l1,
    c_l2 c_l2,
    c_l3 c_l3,
    c_h1 c_h1,
    c_h2 c_h2,
    c_h3 c_h3,
    fine fine,
    f_l1 f_l1,
    f_l2 f_l2,
    f_l3 f_l3,
    f_h1 f_h1,
    f_h2 f_h2,
    f_h3 f_h3,
    lo lo,
    hi hi,
    lrat lrat
) 
```



### <a href="#function-fort-edge-interp" id="function-fort-edge-interp">function fort\_edge\_interp </a>


```cpp
subroutine fort_edge_interp (
    flo flo,
    fhi fhi,
    nc nc,
    ratio ratio,
    dir dir,
    fine fine,
    fine_l0 fine_l0,
    fine_l1 fine_l1,
    fine_l2 fine_l2,
    fine_h0 fine_h0,
    fine_h1 fine_h1,
    fine_h2 fine_h2
) 
```



### <a href="#function-fort-pc-edge-interp" id="function-fort-pc-edge-interp">function fort\_pc\_edge\_interp </a>


```cpp
subroutine fort_pc_edge_interp (
    lo lo,
    hi hi,
    nc nc,
    ratio ratio,
    dir dir,
    crse crse,
    crse_l0 crse_l0,
    crse_l1 crse_l1,
    crse_l2 crse_l2,
    crse_h0 crse_h0,
    crse_h1 crse_h1,
    crse_h2 crse_h2,
    fine fine,
    fine_l0 fine_l0,
    fine_l1 fine_l1,
    fine_l2 fine_l2,
    fine_h0 fine_h0,
    fine_h1 fine_h1,
    fine_h2 fine_h2
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/Nyx/axionyx/Source/Gravity/Gravity_3d.f90`