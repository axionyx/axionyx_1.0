
# File update.f90

[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**update.f90**](update_8f90.md)

[Go to the documentation of this file.](update_8f90.md) 


````cpp
subroutine update(dt, e_in, e_out, rpar)  bind(C,name="fort_update_eos")

      use amrex_error_module, only : amrex_abort
      use amrex_fort_module , only : rt => amrex_real
!      use fundamental_constants_module, only: e_to_cgs, density_to_cgs, &
!                                              heat_from_cgs, mp_over_kb
      use eos_module, only: iterate_ne, &
                            nyx_eos_t_given_re, nyx_eos_given_rt
      use atomic_rates_module, ONLY: tcoolmin, tcoolmax, ncooltab, deltat, &
                                     mproton, xhydrogen, &
                                     uvb_density_a, uvb_density_b, mean_rhob, &
                                     betah0, betahe0, betahep, betaff1, betaff4, &
                                     rechp, rechep, rechepp, &
                                     eh0, ehe0, ehep,&
                                     this_z, yhelium, boltzmann, mproton
      use meth_params_module, only: gamma_minus_1
      use vode_aux_module       , only: z_vode,&! rho_vode, T_vode, ne_vode, &                                                                                                                                     
           jh_vode, jhe_vode, i_vode, j_vode, k_vode, fn_vode, nr_vode
      use reion_aux_module   , only: zhi_flash, zheii_flash, flash_h, flash_he, &
           t_zhi, t_zheii, inhomogeneous_on
      use fundamental_constants_module
      use comoving_module, only: comoving_h, comoving_omb
      use comoving_nd_module, only: fort_integrate_comoving_a
      use amrex_constants_module, only: m_pi
      
      use, intrinsic :: iso_c_binding
      
      real(c_double), value, intent(in) :: dt
      real(c_double), intent(inout) :: e_in(1)
      real(c_double), intent(inout) :: rpar(4)
      real(c_double), intent(  out) :: e_out(1)

      real(rt), parameter :: compt_c = 1.01765467d-37, t_cmb = 2.725d0

      real(rt) :: logT, tmp, fhi, flo
      real(rt) :: ahp, ahep, ahepp, ad, geh0, gehe0, gehep
      real(rt) :: bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp
      real(rt) :: lambda_c, lambda_ff, lambda, heat
      real(rt) :: rho, U, a, rho_heat
      real(rt) :: nh, nh0, nhp, nhe0, nhep, nhepp
      real(rt) :: mu, c_ne, c_lambda, lambda_ref, alpha_ref, t_in
      integer :: j
      integer :: print_radius
      CHARACTER(LEN=80) :: FMT

      real(rt) :: z, z_end, a_end, H_reion_z, He_reion_z
      real(rt) :: T_H, T_He
      real(rt) :: species(5)
!      real(rt) :: z_vode, rho_vode, T_vode, ne_vode                                                                                                                                                                
      real(rt) :: rho_vode, T_vode, ne_vode

      t_vode=rpar(1)
      ne_vode=rpar(2)
      rho_vode=rpar(3)
      z_vode=rpar(4)                                                                                                                                                                                               

      fn_vode=fn_vode+1;

      z = 1.d0/a - 1.d0
      call fort_integrate_comoving_a(a, a_end, dt)
      z_end = 1.0d0/a_end - 1.0d0

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

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT

      if (e_in(1) .lt. 0.d0) &
         e_in(1) = tiny(e_in(1))

     ! Converts from code units to CGS                                                                                                                                                                              
      rho = rho_vode * density_to_cgs * (1.0d0+z_vode)**3
      u = e_in(1) * e_to_cgs
      nh  = rho*xhydrogen/mproton

      call iterate_ne(jh_vode, jhe_vode, z_vode, u, t_vode, nh, ne_vode, nh0, nhp, nhe0, nhep, nhepp)

      ne_vode = nh * ne_vode
      nh0   = nh * nh0
      nhp   = nh * nhp
      nhe0  = nh * nhe0
      nhep  = nh * nhep
      nhepp = nh * nhepp

      logt = dlog10(t_vode)

      if( logt .lt.4 .or. logt .gt. 8) then
 !        print*,logT
         call vode_wrapper(dt,rho_vode,t_vode,ne_vode,e_in, &
              t_vode ,ne_vode ,e_out)

         rpar(1)=t_vode
         rpar(2)=ne_vode
         rpar(3)=rho_vode
         rpar(4)=z_vode

      if (e_out(1) .lt. 0.d0) then
         !$OMP CRITICAL
         print *,'negative e exiting strang integration ',z, rho/mean_rhob, t_vode, e_in, e_out
         call flush(6)
         !$OMP END CRITICAL
         t_vode  = 10.0
         ne_vode = 0.0
         mu     = (1.0d0+4.0d0*yhelium) / (1.0d0+yhelium+ne_vode)
         e_out  = t_vode / (gamma_minus_1 * mp_over_kb * mu)
!         e_out=e_in
         !                    call amrex_abort('bad e out of strang')
      end if
      else
      ! Temperature floor
      if (logt .ge. tcoolmax) then ! Only free-free and Compton cooling are relevant                                                                                                                                
         lambda_ff = 1.42d-27 * dsqrt(t_vode) * (1.1d0 + 0.34d0*dexp(-(5.5d0 - logt)**2 / 3.0d0)) &
                              * (nhp + 4.0d0*nhepp)*ne_vode
         lambda_c  = compt_c*t_cmb**4 * ne_vode * (t_vode - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4

         mu = (1.0d0+4.0d0*yhelium) / (1.0d0+yhelium+ne_vode)
         t_in  = gamma_minus_1*mproton/boltzmann * e_in(1) * e_to_cgs * mu
         c_ne = 1/(gamma_minus_1*mproton/boltzmann * e_to_cgs * mu)
         a = 1.d0 / (1.d0 + z_vode)
         c_lambda = heat_from_cgs/(1.0d0+z_vode)**4 / rho_vode / a
!!         energy  = (-lambda_ff -lambda_c) * heat_from_cgs/(1.0d0+z_vode)**4

         ! Convert to the actual term to be used in e_out = e_in + dt*energy                                                                                                                                        
!!         energy  = energy / rho_vode * (1.0d0+z_vode)
         ne_vode = ne_vode / nh

         lambda_ref = 1.42d-27
         lambda_ref = 3.8e-26
         alpha_ref = .34d0!4/3.d0

         lambda_ref = 1.33d-19
         alpha_ref = -.5
         
         e_out = (- lambda_c)*heat_from_cgs/(1.0d0+z_vode)**4 * dt
         e_out = e_out + c_ne*(dt*(c_lambda/c_ne*lambda_ref*(-alpha_ref+1))+t_in**(alpha_ref+1))**(1/(-alpha_ref+1))

!         print*, 'long',e_out , c_ne,dt,c_lambda,c_ne,lambda_ref,(-alpha_ref+1),t_in**(alpha_ref+1)
         ! Convert to the actual term to be used in e_out = e_in + dt*energy
         a = 1.d0 / (1.d0 + z_vode)
         e_out = e_out / rho_vode / a
!       print *, 'enr = ', energy, 'at (i,j,k) ',i_vode,j_vode,k_vode                                                                                                                                               
!       print *, 'rho_heat = ', rho_heat, 'at (i,j,k) ',i_vode,j_vode,k_vode                                                                                                                                        
!      print(FMT), 'cfrh:',fn_vode,e_in,rho_vode,T_vode,rpar(1)                                                                                                                                                     
!       print *, 'rho = ', rho_vode, 'at (i,j,k) ',i_vode,j_vode,k_vode                                                                                                                                             
         rpar(1)=t_vode
         rpar(2)=ne_vode
         rpar(3)=rho_vode
         rpar(4)=z_vode
      if (e_out(1) .lt. 0.d0) then
         !$OMP CRITICAL
         print *,'negative e exiting strang integration ',z, rho/mean_rhob, t_vode, e_in, e_out
         call flush(6)
         !$OMP END CRITICAL
         t_vode  = 10.0
         ne_vode = 0.0
         mu     = (1.0d0+4.0d0*yhelium) / (1.0d0+yhelium+ne_vode)
         e_out  = t_vode / (gamma_minus_1 * mp_over_kb * mu)
         e_out=e_in
         !                    call amrex_abort('bad e out of strang')
      end if
!      e_out=1.d50

      z=z_vode
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
         
         return
      end if
      
      if (logt .le. tcoolmin)  logt = tcoolmin + 0.5d0*deltat

      lambda_ref = 1.33d-19
      alpha_ref = -.5
      
      lambda_c = compt_c*t_cmb**4*ne_vode*(t_vode - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4   ! Compton cooling                                                                                                   
      lambda = lambda + lambda_c

      ! Heating terms                                                                                                                                                                                               
      heat = jh_vode*nh0*eh0 + jh_vode*nhe0*ehe0 + jhe_vode*nhep*ehep
      rho_heat = uvb_density_a * (rho_vode/mean_rhob)**uvb_density_b
      heat = rho_heat*heat

      ! Convert back to code units

      mu = (1.0d0+4.0d0*yhelium) / (1.0d0+yhelium+ne_vode)
      t_in  = gamma_minus_1*mproton/boltzmann * e_in(1) * e_to_cgs * mu
      c_ne = 1/(gamma_minus_1*mproton/boltzmann * e_to_cgs * mu)
      a = 1.d0 / (1.d0 + z_vode)
      c_lambda = heat_from_cgs/(1.0d0+z_vode)**4 / rho_vode / a
      
      ne_vode     = ne_vode / nh
      e_out = (heat - lambda_c)*heat_from_cgs/(1.0d0+z_vode)**4 * dt
      e_out = e_out + c_ne*(dt*(c_lambda/c_ne*lambda_ref*(-alpha_ref+1))+t_in**(alpha_ref+1))**(1/(-alpha_ref+1))
      ! Convert to the actual term to be used in e_out = e_in + dt*energy                                                                                                                                           
      e_out = e_out / rho_vode / a
!      print(FMT), 'dfrh:',fn_vode,e_in,rho_vode,T_vode,rpar(1)                                                                                                                                                     
      rpar(1)=t_vode
      rpar(2)=ne_vode
      rpar(3)=rho_vode
      rpar(4)=z_vode

      if (e_out(1) .lt. 0.d0) then
         !$OMP CRITICAL
         print *,'negative e exiting strang integration ',z, rho/mean_rhob, t_vode, e_in, e_out
         call flush(6)
         !$OMP END CRITICAL
         t_vode  = 10.0
         ne_vode = 0.0
         mu     = (1.0d0+4.0d0*yhelium) / (1.0d0+yhelium+ne_vode)
         e_out  = t_vode / (gamma_minus_1 * mp_over_kb * mu)
         e_out=e_in
         !                    call amrex_abort('bad e out of strang')
      end if
      endif

      a = 1.d0 / (1.d0 + z_vode)
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
      
end subroutine update
````

