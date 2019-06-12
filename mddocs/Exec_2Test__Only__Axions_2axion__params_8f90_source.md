
# File axion\_params.f90

[**File List**](files.md) **>** [**Exec**](dir_43a12cefb7942b6f49b5b628aafd3192.md) **>** [**Test\_Only\_Axions**](dir_eb24725df855cf6c732a19e4912f662a.md) **>** [**axion\_params.f90**](Exec_2Test__Only__Axions_2axion__params_8f90.md)

[Go to the documentation of this file.](Exec_2Test__Only__Axions_2axion__params_8f90.md) 


````cpp

module axion_params_module

  ! These are only used for the axion model.
  ! They have to be set/overwritten in Prob_3d.f90.
  
  double precision, save :: m_tt = 2.5d0 !particle mass in units of 10^(-22) eV
  double precision, save :: meandens !background axion density, set in ca_initdata

  double complex, parameter :: ii = (0., 1.)   

  double precision, save :: epsilon_L = 0.3d0  ! needed for vel. tagging (cf. Loehner 1987)
  double precision, save :: mindens   = 1.0d+3  ! velocity tagging only in regions where (density > mindens)
  double precision, save :: critvalue = 0.3d0  ! tagging regions where error indicator 0<=err_ind<1 (cf. Loehner 1987 eq.4) is larger then critvalue

end module axion_params_module


!! for now this is one is called from Prob_3d.f90
!subroutine fort_set_nu_small_dens
!  
!  use comoving_module, only : comoving_h, comoving_OmNu
!  use fundamental_constants_module, only : Gconst
!  use bl_constants_module, only : M_PI
!  use nu_params_module, only : nu_small_dens
!
!  implicit none
!
!  double precision :: rhoNuBar
!
!  ! this is the mean neutrino density from \Omega_\nu
!  rhoNuBar = 3 * (comoving_h*100) * (comoving_h*100) * comoving_OmNu &
!             / 8.d0 / M_PI / Gconst
!
!  nu_small_dens = 1d-6*rhoNuBar
!
!end subroutine fort_set_nu_small_dens




````

