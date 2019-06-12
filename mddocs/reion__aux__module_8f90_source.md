
# File reion\_aux\_module.f90

[**File List**](files.md) **>** [**EOS**](dir_2a6406f09975eea078703cc63b0e3416.md) **>** [**reion\_aux\_module.f90**](reion__aux__module_8f90.md)

[Go to the documentation of this file.](reion__aux__module_8f90.md) 


````cpp
module reion_aux_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! Global variables (re)set on inputs
  real(rt), save :: zhi_flash=-1.0, zheii_flash=-1.0, t_zhi=0.0, t_zheii=0.0
  logical, save  :: flash_h=.false., flash_he=.false., inhomogeneous_on=.false.

end module reion_aux_module
````

