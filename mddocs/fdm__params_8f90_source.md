
# File fdm\_params.f90

[**File List**](files.md) **>** [**FDM**](dir_43b815edcf2a06ee60d8a45cc6c77fb8.md) **>** [**fdm\_params.f90**](fdm__params_8f90.md)

[Go to the documentation of this file.](fdm__params_8f90.md) 


````cpp
module fdm_params_module

  use amrex_fort_module, only : rt => amrex_real
  
  real(rt), save :: m_tt
  real(rt), save :: hbaroverm
  real(rt), save :: theta_fdm
  real(rt), save :: sigma_fdm
  real(rt), save :: gamma_fdm
  real(rt), save :: meandens !background fdm density, set in ca_initdata
  real(rt), save :: a
  ! integer, save :: wkb_approx 

  double complex, parameter :: ii = (0., 1.)   

  real(rt), save :: epsilon_L = 0.3d0  ! needed for vel. tagging (cf. Loehner 1987)
  real(rt), save :: mindens   = 1.0d+3  ! velocity tagging only in regions where (density > mindens)
  real(rt), save :: critvalue = 0.3d0  ! tagging regions where error indicator 0<=err_ind<1 (cf. Loehner 1987 eq.4) is larger then critvalue

end module fdm_params_module
````

