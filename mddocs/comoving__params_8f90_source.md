
# File comoving\_params.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**comoving\_params.f90**](comoving__params_8f90.md)

[Go to the documentation of this file.](comoving__params_8f90.md) 


````cpp
module comoving_module

    use amrex_fort_module, only : rt => amrex_real

    integer,  save :: comoving_type
    real(rt), save :: comoving_OmM, comoving_OmB, comoving_OmN, comoving_h, comoving_OmL
    !for axions/scalar field dark matter:
    real(rt), save :: comoving_OmAx

end module comoving_module
````

