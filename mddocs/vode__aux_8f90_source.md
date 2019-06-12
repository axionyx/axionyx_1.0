
# File vode\_aux.f90

[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**vode\_aux.f90**](vode__aux_8f90.md)

[Go to the documentation of this file.](vode__aux_8f90.md) 


````cpp

! This module stores the extra parameters for the VODE calls.

module vode_aux_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt), save :: z_vode
  real(rt), save :: rho_vode, T_vode, ne_vode, rho_init_vode, rho_src_vode, rhoe_src_vode, e_src_vode
  real(rt), dimension(:), allocatable, save :: rho_vode_vec, T_vode_vec, ne_vode_vec
  integer , save :: JH_vode, JHe_vode, i_vode, i_point, j_point, k_point, j_vode, k_vode, fn_vode, NR_vode
  logical,  save :: firstcall
  !$OMP THREADPRIVATE (rho_vode, rho_vode_vec, T_vode, T_vode_vec, ne_vode, ne_vode_vec, JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, firstcall, rho_init_vode, rho_src_vode,rhoe_src_vode, e_src_vode, i_point, j_point, k_point)

end module vode_aux_module
````

