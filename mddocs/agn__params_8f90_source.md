
# File agn\_params.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**agn\_params.f90**](agn__params_8f90.md)

[Go to the documentation of this file.](agn__params_8f90.md) 


````cpp
module agn_params_module

    use amrex_fort_module, only : rt => amrex_real

    real(rt), save :: l_merge
    logical, save :: cutoff_vel
    real(rt), save :: eps_rad, eps_coupling, T_min, bondi_boost, max_frac_removed, frac_kinetic, eps_kinetic

end module agn_params_module
````

