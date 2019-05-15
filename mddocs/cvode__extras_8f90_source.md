
# File cvode\_extras.f90

[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**cvode\_extras.f90**](cvode__extras_8f90.md)

[Go to the documentation of this file.](cvode__extras_8f90.md) 


````cpp
module cvode_extras
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  contains

    subroutine ode_eos_setup(a,half_dt) &
      bind(C,name="fort_ode_eos_setup")
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : nvar, urho, ueden, ueint, &
                                   ndiag, temp_comp, ne_comp, zhi_comp, &
                                   gamma_minus_1
    use amrex_constants_module, only: m_pi
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_t_given_re, nyx_eos_given_rt
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_omb
    use comoving_nd_module, only: fort_integrate_comoving_a
    use atomic_rates_module, only: yhelium
    use vode_aux_module    , only: jh_vode, jhe_vode, z_vode, i_vode, j_vode, k_vode
    use reion_aux_module   , only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                   t_zhi, t_zheii, inhomogeneous_on
    real(rt), intent(inout) :: a
    real(rt), intent(inout) :: half_dt
    real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
     real(rt) :: mu, mean_rhob, T_H, T_He
    real(rt) :: species(5)
    real(rt) :: rho_vode, T_vode, ne_vode
    z = 1.d0/a - 1.d0
    call fort_integrate_comoving_a(a, a_end, half_dt)
    z_end = 1.0d0/a_end - 1.0d0
    print*, z,z_end,a,a_end,half_dt

    mean_rhob = comoving_omb * 3.d0*(comoving_h*100.d0)**2 / (8.d0*m_pi*gconst)

    ! Flash reionization?
    if ((flash_h .eqv. .true.) .and. (z .gt. zhi_flash)) then
       jh_vode = 0
    else
       jh_vode = 1
    endif
    if ((flash_he .eqv. .true.) .and. (z .gt. zheii_flash)) then
       jhe_vode = 0
    else
       jhe_vode = 1
    endif

    if (flash_h ) h_reion_z  = zhi_flash
    if (flash_he) he_reion_z = zheii_flash

    if (inhomogeneous_on) then
       stop "Do not currently support inhomogenous_on with box"
       !H_reion_z = diag_eos(i,j,k,ZHI_COMP)
       if (z .gt. h_reion_z) then
          jh_vode = 0
       else
          jh_vode = 1
       endif
    endif

    end subroutine ode_eos_setup

    subroutine ode_eos_finalize(e_out, rpar, num_eq) &
      bind(C,name="fort_ode_eos_finalize")
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : nvar, urho, ueden, ueint, &
                                   ndiag, temp_comp, ne_comp, zhi_comp, &
                                   gamma_minus_1
    use amrex_constants_module, only: m_pi
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_t_given_re, nyx_eos_given_rt
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_omb
    use comoving_nd_module, only: fort_integrate_comoving_a
    use atomic_rates_module, only: yhelium
    use vode_aux_module    , only: jh_vode, jhe_vode, z_vode, i_vode, j_vode, k_vode
    use reion_aux_module   , only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                   t_zhi, t_zheii, inhomogeneous_on
    integer, value :: num_eq
   real(rt), intent(inout) :: e_out(num_eq)
   real(rt), intent(inout) :: rpar(4*num_eq)

    real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
     real(rt) :: mu, mean_rhob, T_H, T_He
    real(rt) :: species(5)
    real(rt) :: rho_vode, T_vode, ne_vode
    real(rt) :: a

      t_vode=rpar(1)
      ne_vode=rpar(2)
      rho_vode=rpar(3)
      a=1/(rpar(4)+1)


      if (e_out(1) .lt. 0.d0) then
         !$OMP CRITICAL
         print *,'negative e exiting strang integration ',z, rho/mean_rhob, e_out
         call flush(6)
         !$OMP END CRITICAL
         t_vode  = 10.0
         ne_vode = 0.0
         mu     = (1.0d0+4.0d0*yhelium) / (1.0d0+yhelium+ne_vode)
         e_out  = t_vode / (gamma_minus_1 * mp_over_kb * mu)
         !                    call amrex_abort('bad e out of strang')
      end if
      
      ! Update T and ne (do not use stuff computed in f_rhs, per vode manual)
      call nyx_eos_t_given_re(jh_vode, jhe_vode, t_vode, ne_vode, rho_vode, e_out(1), a, species)
      
      ! Instanteneous heating from reionization:
      t_h = 0.0d0
      if (inhomogeneous_on .or. flash_h) then
         if ((h_reion_z  .lt. z) .and. (h_reion_z  .ge. z_end)) t_h  = (1.0d0 - species(2))*max((t_zhi-t_vode), 0.0d0)
      endif
      
      t_he = 0.0d0
      if (flash_he) then
         if ((he_reion_z .lt. z) .and. (he_reion_z .ge. z_end)) t_he = (1.0d0 - species(5))*max((t_zheii-t_vode), 0.0d0)
      endif
      
      if ((t_h .gt. 0.0d0) .or. (t_he .gt. 0.0d0)) then
         t_vode = t_vode + t_h + t_he                            ! For simplicity, we assume
         ne_vode = 1.0d0 + yhelium                              !    completely ionized medium at
         if (t_he .gt. 0.0d0) ne_vode = ne_vode + yhelium        !    this point.  It's a very minor
         mu = (1.0d0+4.0d0*yhelium) / (1.0d0+yhelium+ne_vode)   !    detail compared to the overall approximation.
         e_out  = t_vode / (gamma_minus_1 * mp_over_kb * mu)
         call nyx_eos_t_given_re(jh_vode, jhe_vode, t_vode, ne_vode, rho_vode, e_out(1), a, species)
      endif
      rpar(1)=t_vode
      rpar(2)=ne_vode
      rpar(3)=rho_vode
      rpar(4)=z_vode      

    end subroutine ode_eos_finalize

    integer(c_int) function rhsfnreal(tn, yvec, fvec, rpar, neq) &
           result(ierr) bind(c,name='RhsFnReal')

      use, intrinsic :: iso_c_binding
      implicit none


      real(c_double), value :: tn
      integer(c_int), value :: neq
!      type(c_ptr), value    :: sunvec_y                                                                                                                                                                            
!      type(c_ptr), value    :: sunvec_f                                                                                                                                                                            
!      type(c_ptr), value    :: user_data                                                                                                                                                                           

      ! pointers to data in SUNDAILS vectors                                                                                                                                                                        
      real(c_double) :: yvec(neq)
      real(c_double) :: fvec(neq)
      real(c_double), intent(inout) :: rpar(neq*4)
      real(c_double) :: energy(neq)

!      print*, "r1", rpar(1)                                                                                                                                                                                        
!      print*, "r2", rpar(2)                                                                                                                                                                                        
!      print*, rpar(3)                                                                                                                                                                                              
!      print*, rpar(4)                                                                                                                                                                                              
      call f_rhs_rpar(neq, tn, yvec, fvec, rpar, 0)
!      print*, "after r1", rpar(1)                                                                                                                                                                                  
!      print*, "after r2", rpar(2)                                                                                                                                                                                  
!      print*, "after r3", rpar(3)                                                                                                                                                                                  
!      print*, "after r4", rpar(4)                                                                                                                                                                                  

      ierr = 0
    end function rhsfnreal

end module cvode_extras
````

