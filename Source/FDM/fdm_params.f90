module fdm_params_module

  use amrex_fort_module, only : rt => amrex_real
  
  real(rt), save :: m_tt
  real(rt), save :: hbaroverm
  real(rt), save :: theta_fdm
  real(rt), save :: sigma_fdm
  real(rt), save :: gamma_fdm
  real(rt), save :: meandens !background fdm density, set in ca_initdata
  real(rt), save :: a
  real(rt), save :: ratio_fdm
  ! integer, save :: wkb_approx 

  double complex, parameter :: ii = (0., 1.)   

  real(rt), save :: epsilon_L = 0.3d0  ! needed for vel. tagging (cf. Loehner 1987)
  real(rt), save :: mindens   = 1.0d11  ! velocity tagging only in regions where (density > mindens)
  real(rt), save :: critvalue = 0.3d0  ! tagging regions where error indicator 0<=err_ind<1 (cf. Loehner 1987 eq.4) is larger then critvalue

end module fdm_params_module
