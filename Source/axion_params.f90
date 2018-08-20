
module axion_params_module

  ! These are only used for the neutrino fluid model.
  ! They have to be set/overwritten in Prob_3d.f90.
!  double precision, save :: VISCOSITY      = 0.d0/1.d0
!  double precision, save :: SOUNDSPEED     = 1.d0
!  double precision, save :: SOUNDSPEED_OLD = 1.d0
!  double precision, save :: SOUNDSPEED_NEW = 1.d0
!
!  double precision, save :: nu_small_dens  = -1d-200

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




