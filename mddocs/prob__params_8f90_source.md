
# File prob\_params.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**prob\_params.f90**](prob__params_8f90.md)

[Go to the documentation of this file.](prob__params_8f90.md) 


````cpp

! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_params_module

  implicit none
  integer         , save :: physbc_lo(3)
  integer         , save :: physbc_hi(3)
  integer         , save :: Outflow, Symmetry
  integer         , save :: coord_type

end module prob_params_module
````

