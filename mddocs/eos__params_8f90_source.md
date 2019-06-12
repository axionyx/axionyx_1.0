
# File eos\_params.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**eos\_params.f90**](eos__params_8f90.md)

[Go to the documentation of this file.](eos__params_8f90.md) 


````cpp

! This module stores the runtime EOS species IF they are defined to be constants.  
! These parameter are initialized in set_eos_params().

module eos_params_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), save ::  h_species
  real(rt), save :: he_species

end module eos_params_module
````

