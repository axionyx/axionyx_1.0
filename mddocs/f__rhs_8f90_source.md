
# File f\_rhs.f90

[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**f\_rhs.f90**](f__rhs_8f90.md)

[Go to the documentation of this file.](f__rhs_8f90.md) 


````cpp
subroutine f_rhs_rpar(num_eq, time, e_in, energy, rpar, ipar)

      use amrex_error_module, only : amrex_abort
      use amrex_fort_module , only : rt => amrex_real
      use fundamental_constants_module, only: e_to_cgs, density_to_cgs, &
                                              heat_from_cgs
      use eos_module, only: iterate_ne
      use atomic_rates_module, ONLY: tcoolmin, tcoolmax, ncooltab, deltat, &
                                     mproton, xhydrogen, &
                                     uvb_density_a, uvb_density_b, mean_rhob, &
                                     betah0, betahe0, betahep, betaff1, betaff4, &
                                     rechp, rechep, rechepp, &
                                     eh0, ehe0, ehep

      use vode_aux_module       , only: z_vode,&! rho_vode, T_vode, ne_vode, &                                                                                                                                     
                                        jh_vode, jhe_vode, i_vode, j_vode, k_vode, fn_vode, nr_vode

      integer, intent(in)             :: num_eq, ipar
      real(rt), intent(inout) :: e_in(num_eq)
      real(rt), intent(in   ) :: time
      real(rt), intent(inout) :: rpar(4*num_eq)
      real(rt), intent(  out) :: energy(num_eq)

      real(rt), parameter :: compt_c = 1.01765467d-37, t_cmb = 2.725d0

      real(rt) :: logT, tmp, fhi, flo
      real(rt) :: ahp, ahep, ahepp, ad, geh0, gehe0, gehep
      real(rt) :: bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp
      real(rt) :: lambda_c, lambda_ff, lambda, heat
      real(rt) :: rho, U, a, rho_heat
      real(rt) :: nh, nh0, nhp, nhe0, nhep, nhepp
      integer :: j
      integer :: print_radius
      CHARACTER(LEN=80) :: FMT
!      real(rt) :: z_vode, rho_vode, T_vode, ne_vode                                                                                                                                                                
      real(rt) :: rho_vode, T_vode, ne_vode

      t_vode=rpar(1)
      ne_vode=rpar(2)
      rho_vode=rpar(3)
!      z_vode=rpar(4)                                                                                                                                                                                               

      fn_vode=fn_vode+1;

!      FMT = "(A6, I4, ES15.5, ES15.5E3, ES15.5, ES15.5)"                                                                                                                                                           
!      print(FMT), 'frh:',fn_vode,e_in,rho_vode,T_vode,rpar(1)                                                                                                                                                      

      if (e_in(1) .lt. 0.d0) &
         e_in(1) = tiny(e_in(1))

     ! Converts from code units to CGS                                                                                                                                                                              
      rho = rho_vode * density_to_cgs * (1.0d0+z_vode)**3
      u = e_in(1) * e_to_cgs
      nh  = rho*xhydrogen/mproton
!      print*, "rho = ", rho                                                                                                                                                                                        
!      print*, "nh = ", nh                                                                                                                                                                                          
!      nh = 1.85e-2                                                                                                                                                                                                 
!      print*, "Replacing with nh = ",nh                                                                                                                                                                            
!      print*, "XHYDROGEN = ", XHYDROGEN                                                                                                                                                                            
!      print*, "MPROTON = ", MPROTON                                                                                                                                                                                
!      print*, "density_to_cgs = ", density_to_cgs                                                                                                                                                                  
!      print*, "e_to_cgs = ", e_to_cgs                                                                                                                                                                              

      if (time .gt. 1) then
         print *,'TIME INTO F_RHS ',time
         print *,'AT              ',i_vode,j_vode,k_vode
         call amrex_abort("TOO BIG TIME IN F_RHS") 
      end if
!      print(FMT), 'afrh:',fn_vode,e_in,rho_vode,T_vode,rpar(1)                                                                                                                                                     
      ! Get gas temperature and individual ionization species                                                                                                                                                       
      ! testing different memory structures                                                                                                                                                                         
!      NR_vode=0                                                                                                                                                                                                    
      call iterate_ne(jh_vode, jhe_vode, z_vode, u, t_vode, nh, ne_vode, nh0, nhp, nhe0, nhep, nhepp)
!      print(FMT), 'bfrh:',fn_vode,e_in,rho_vode,T_vode,time                                                                                                                                                        
      ! Convert species to CGS units:                                                                                                                                                                               
      ne_vode = nh * ne_vode
      nh0   = nh * nh0
      nhp   = nh * nhp
      nhe0  = nh * nhe0
      nhep  = nh * nhep
      nhepp = nh * nhepp

      logt = dlog10(t_vode)
      if (logt .ge. tcoolmax) then ! Only free-free and Compton cooling are relevant                                                                                                                                
         lambda_ff = 1.42d-27 * dsqrt(t_vode) * (1.1d0 + 0.34d0*dexp(-(5.5d0 - logt)**2 / 3.0d0)) &
                              * (nhp + 4.0d0*nhepp)*ne_vode
         lambda_c  = compt_c*t_cmb**4 * ne_vode * (t_vode - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4

         energy  = (-lambda_ff -lambda_c) * heat_from_cgs/(1.0d0+z_vode)**4

         ! Convert to the actual term to be used in e_out = e_in + dt*energy                                                                                                                                        
         energy  = energy / rho_vode * (1.0d0+z_vode)
         ne_vode = ne_vode / nh

!       print *, 'enr = ', energy, 'at (i,j,k) ',i_vode,j_vode,k_vode                                                                                                                                               
!       print *, 'rho_heat = ', rho_heat, 'at (i,j,k) ',i_vode,j_vode,k_vode                                                                                                                                        
!      print(FMT), 'cfrh:',fn_vode,e_in,rho_vode,T_vode,rpar(1)                                                                                                                                                     
!       print *, 'rho = ', rho_vode, 'at (i,j,k) ',i_vode,j_vode,k_vode                                                                                                                                             
         rpar(1)=t_vode
         rpar(2)=ne_vode
         rpar(3)=rho_vode
         rpar(4)=z_vode
         return
      end if

      ! Temperature floor                                                                                                                                                                                           
      if (logt .le. tcoolmin)  logt = tcoolmin + 0.5d0*deltat

      ! Interpolate rates                                                                                                                                                                                           
      tmp = (logt-tcoolmin)/deltat
      j = int(tmp)
      fhi = tmp - j
      flo = 1.0d0 - fhi
      j = j + 1 ! F90 arrays start with 1              
      bh0   = flo*betah0(j) + fhi*betah0(j+1)
      bhe0  = flo*betahe0(j) + fhi*betahe0(j+1)
      bhep  = flo*betahep(j) + fhi*betahep(j+1)
      bff1  = flo*betaff1(j) + fhi*betaff1(j+1)
      bff4  = flo*betaff4(j) + fhi*betaff4(j+1)
      rhp   = flo*rechp(j) + fhi*rechp(j+1)
      rhep  = flo*rechep(j) + fhi*rechep(j+1)
      rhepp = flo*rechepp(j) + fhi*rechepp(j+1)

      ! Cooling:                                                                                                                                                                                                    
      lambda = ( bh0*nh0 + bhe0*nhe0 + bhep*nhep + &
                 rhp*nhp + rhep*nhep + rhepp*nhepp + &
                 bff1*(nhp+nhep) + bff4*nhepp ) * ne_vode

      lambda_c = compt_c*t_cmb**4*ne_vode*(t_vode - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4   ! Compton cooling                                                                                                   
      lambda = lambda + lambda_c

      ! Heating terms                                                                                                                                                                                               
      heat = jh_vode*nh0*eh0 + jh_vode*nhe0*ehe0 + jhe_vode*nhep*ehep
      rho_heat = uvb_density_a * (rho_vode/mean_rhob)**uvb_density_b
      heat = rho_heat*heat

      ! Convert back to code units                                                                                                                                                                                  
      ne_vode     = ne_vode / nh
      energy = (heat - lambda)*heat_from_cgs/(1.0d0+z_vode)**4
      ! Convert to the actual term to be used in e_out = e_in + dt*energy                                                                                                                                           
      a = 1.d0 / (1.d0 + z_vode)
      energy = energy / rho_vode / a
!      print(FMT), 'dfrh:',fn_vode,e_in,rho_vode,T_vode,rpar(1)                                                                                                                                                     
      rpar(1)=t_vode
      rpar(2)=ne_vode
      rpar(3)=rho_vode
      rpar(4)=z_vode
!      print(FMT), 'efrh:',fn_vode,e_in,rho_vode,T_vode,rpar(1)                                                                                                                                                     
!       print *, 'energy = ', energy, 'at (i,j,k) ',i_vode,j_vode,k_vode                                                                                                                                            
!       print *, 'rho_heat = ', rho_heat, 'at (i,j,k) ',i_vode,j_vode,k_vode                                                                                                                                        
!       print *, 'rho_vd = ', rho_vode, 'at (i,j,k) ',i_vode,j_vode,k_vode         
end subroutine f_rhs_rpar

subroutine f_rhs_split(num_eq, time, y_in, yp_out, rpar, ipar)

      use amrex_error_module, only : amrex_abort
      use amrex_fort_module, only : rt => amrex_real
      use fundamental_constants_module, only: e_to_cgs, density_to_cgs, & 
                                              heat_from_cgs
      use eos_module, only: iterate_ne
      use atomic_rates_module, ONLY: tcoolmin, tcoolmax, ncooltab, deltat, &
                                     mproton, xhydrogen, &
                                     uvb_density_a, uvb_density_b, mean_rhob, &
                                     betah0, betahe0, betahep, betaff1, betaff4, &
                                     rechp, rechep, rechepp, &
                                     eh0, ehe0, ehep

      use vode_aux_module       , only: z_vode, rho_vode, t_vode, ne_vode, &
                                        jh_vode, jhe_vode, i_vode, j_vode, k_vode, fn_vode, nr_vode, &
                                        rho_init_vode, e_src_vode, rho_src_vode

      integer, intent(in)             :: num_eq, ipar
      real(rt), intent(inout) :: y_in(num_eq)
      real(rt), intent(in   ) :: time
      real(rt), intent(in   ) :: rpar
      real(rt), intent(  out) :: yp_out(num_eq)

      real(rt), parameter :: compt_c = 1.01765467d-37, t_cmb = 2.725d0


      real(rt) :: e_in(1)
      real(rt) :: energy
      real(rt) :: rho_in
      real(rt) :: logT, tmp, fhi, flo
      real(rt) :: ahp, ahep, ahepp, ad, geh0, gehe0, gehep
      real(rt) :: bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp
      real(rt) :: lambda_c, lambda_ff, lambda, heat
      real(rt) :: rho, U, a, rho_heat
      real(rt) :: nh, nh0, nhp, nhe0, nhep, nhepp
      integer :: j

      fn_vode=fn_vode+1;

      e_in = y_in(1)
      rho_vode = y_in(2)
!      rho_vode = rho_init_vode + time * rho_src_vode

      if (e_in(1) .lt. 0.d0) &
         e_in(1) = tiny(e_in(1))

     ! Converts from code units to CGS
      rho = rho_vode * density_to_cgs * (1.0d0+z_vode)**3
        u = e_in(1) * e_to_cgs
      nh  = rho*xhydrogen/mproton

      if (time .gt. 1) then
         print *,'TIME INTO F_RHS ',time
         print *,'AT              ',i_vode,j_vode,k_vode
         call amrex_abort("TOO BIG TIME IN F_RHS")
      end if

      ! Get gas temperature and individual ionization species
      ! testing different memory structures
!      NR_vode=0
      call iterate_ne(jh_vode, jhe_vode, z_vode, u, t_vode, nh, ne_vode, nh0, nhp, nhe0, nhep, nhepp)

      ! Convert species to CGS units: 
      ne_vode = nh * ne_vode
      nh0   = nh * nh0
      nhp   = nh * nhp
      nhe0  = nh * nhe0
      nhep  = nh * nhep
      nhepp = nh * nhepp

      logt = dlog10(t_vode)
      if (logt .ge. tcoolmax) then ! Only free-free and Compton cooling are relevant
         lambda_ff = 1.42d-27 * dsqrt(t_vode) * (1.1d0 + 0.34d0*dexp(-(5.5d0 - logt)**2 / 3.0d0)) &
                              * (nhp + 4.0d0*nhepp)*ne_vode
         lambda_c  = compt_c*t_cmb**4 * ne_vode * (t_vode - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4

         energy  = (-lambda_ff -lambda_c) * heat_from_cgs/(1.0d0+z_vode)**4

         ! Convert to the actual term to be used in e_out = e_in + dt*energy
         energy  = energy / rho_vode * (1.0d0+z_vode)
         ne_vode = ne_vode / nh
         yp_out(1) = energy + e_src_vode
         yp_out(2) = rho_src_vode
         return
      end if

      ! Temperature floor
      if (logt .le. tcoolmin)  logt = tcoolmin + 0.5d0*deltat

      ! Interpolate rates
      tmp = (logt-tcoolmin)/deltat
      j = int(tmp)
      fhi = tmp - j
      flo = 1.0d0 - fhi
      j = j + 1 ! F90 arrays start with 1

      bh0   = flo*betah0(j) + fhi*betah0(j+1)
      bhe0  = flo*betahe0(j) + fhi*betahe0(j+1)
      bhep  = flo*betahep(j) + fhi*betahep(j+1)
      bff1  = flo*betaff1(j) + fhi*betaff1(j+1)
      bff4  = flo*betaff4(j) + fhi*betaff4(j+1)
      rhp   = flo*rechp(j) + fhi*rechp(j+1)
      rhep  = flo*rechep(j) + fhi*rechep(j+1)
      rhepp = flo*rechepp(j) + fhi*rechepp(j+1)

      ! Cooling: 
      lambda = ( bh0*nh0 + bhe0*nhe0 + bhep*nhep + &
                 rhp*nhp + rhep*nhep + rhepp*nhepp + &
                 bff1*(nhp+nhep) + bff4*nhepp ) * ne_vode

      lambda_c = compt_c*t_cmb**4*ne_vode*(t_vode - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4   ! Compton cooling
      lambda = lambda + lambda_c

      ! Heating terms
      heat = jh_vode*nh0*eh0 + jh_vode*nhe0*ehe0 + jhe_vode*nhep*ehep
      rho_heat = uvb_density_a * (rho_vode/mean_rhob)**uvb_density_b
      heat = rho_heat*heat

      ! Convert back to code units
      ne_vode     = ne_vode / nh
      energy = (heat - lambda)*heat_from_cgs/(1.0d0+z_vode)**4

      ! Convert to the actual term to be used in e_out = e_in + dt*energy
      a = 1.d0 / (1.d0 + z_vode)
      energy = (energy) / rho_vode / a

      yp_out(1) = energy + e_src_vode
      yp_out(2) = rho_src_vode


end subroutine f_rhs_split


subroutine f_rhs(num_eq, time, e_in, energy, rpar, ipar)

      use amrex_error_module, only : amrex_abort
      use amrex_fort_module, only : rt => amrex_real
      use fundamental_constants_module, only: e_to_cgs, density_to_cgs, & 
                                              heat_from_cgs
      use eos_module, only: iterate_ne
      use atomic_rates_module, ONLY: tcoolmin, tcoolmax, ncooltab, deltat, &
                                     mproton, xhydrogen, &
                                     uvb_density_a, uvb_density_b, mean_rhob, &
                                     betah0, betahe0, betahep, betaff1, betaff4, &
                                     rechp, rechep, rechepp, &
                                     eh0, ehe0, ehep

      use vode_aux_module       , only: z_vode, rho_vode, t_vode, ne_vode, &
                                        jh_vode, jhe_vode, i_vode, j_vode, k_vode, fn_vode, nr_vode

      integer, intent(in)             :: num_eq, ipar
      real(rt), intent(inout) :: e_in(num_eq)
      real(rt), intent(in   ) :: time
      real(rt), intent(in   ) :: rpar
      real(rt), intent(  out) :: energy

      real(rt), parameter :: compt_c = 1.01765467d-37, t_cmb = 2.725d0

      real(rt) :: logT, tmp, fhi, flo
      real(rt) :: ahp, ahep, ahepp, ad, geh0, gehe0, gehep
      real(rt) :: bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp
      real(rt) :: lambda_c, lambda_ff, lambda, heat
      real(rt) :: rho, U, a, rho_heat
      real(rt) :: nh, nh0, nhp, nhe0, nhep, nhepp
      integer :: j

      fn_vode=fn_vode+1;

      if (e_in(1) .lt. 0.d0) &
         e_in(1) = tiny(e_in(1))

     ! Converts from code units to CGS
      rho = rho_vode * density_to_cgs * (1.0d0+z_vode)**3
        u = e_in(1) * e_to_cgs
      nh  = rho*xhydrogen/mproton

      if (time .gt. 1) then
         print *,'TIME INTO F_RHS ',time
         print *,'AT              ',i_vode,j_vode,k_vode
         call amrex_abort("TOO BIG TIME IN F_RHS")
      end if

      ! Get gas temperature and individual ionization species
      ! testing different memory structures
!      NR_vode=0
      call iterate_ne(jh_vode, jhe_vode, z_vode, u, t_vode, nh, ne_vode, nh0, nhp, nhe0, nhep, nhepp)

      ! Convert species to CGS units: 
      ne_vode = nh * ne_vode
      nh0   = nh * nh0
      nhp   = nh * nhp
      nhe0  = nh * nhe0
      nhep  = nh * nhep
      nhepp = nh * nhepp

      logt = dlog10(t_vode)
      if (logt .ge. tcoolmax) then ! Only free-free and Compton cooling are relevant
         lambda_ff = 1.42d-27 * dsqrt(t_vode) * (1.1d0 + 0.34d0*dexp(-(5.5d0 - logt)**2 / 3.0d0)) &
                              * (nhp + 4.0d0*nhepp)*ne_vode
         lambda_c  = compt_c*t_cmb**4 * ne_vode * (t_vode - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4

         energy  = (-lambda_ff -lambda_c) * heat_from_cgs/(1.0d0+z_vode)**4

         ! Convert to the actual term to be used in e_out = e_in + dt*energy
         energy  = energy / rho_vode * (1.0d0+z_vode)
         ne_vode = ne_vode / nh
         return
      end if

      ! Temperature floor
      if (logt .le. tcoolmin)  logt = tcoolmin + 0.5d0*deltat

      ! Interpolate rates
      tmp = (logt-tcoolmin)/deltat
      j = int(tmp)
      fhi = tmp - j
      flo = 1.0d0 - fhi
      j = j + 1 ! F90 arrays start with 1

      bh0   = flo*betah0(j) + fhi*betah0(j+1)
      bhe0  = flo*betahe0(j) + fhi*betahe0(j+1)
      bhep  = flo*betahep(j) + fhi*betahep(j+1)
      bff1  = flo*betaff1(j) + fhi*betaff1(j+1)
      bff4  = flo*betaff4(j) + fhi*betaff4(j+1)
      rhp   = flo*rechp(j) + fhi*rechp(j+1)
      rhep  = flo*rechep(j) + fhi*rechep(j+1)
      rhepp = flo*rechepp(j) + fhi*rechepp(j+1)

      ! Cooling: 
      lambda = ( bh0*nh0 + bhe0*nhe0 + bhep*nhep + &
                 rhp*nhp + rhep*nhep + rhepp*nhepp + &
                 bff1*(nhp+nhep) + bff4*nhepp ) * ne_vode

      lambda_c = compt_c*t_cmb**4*ne_vode*(t_vode - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4   ! Compton cooling
      lambda = lambda + lambda_c

      ! Heating terms
      heat = jh_vode*nh0*eh0 + jh_vode*nhe0*ehe0 + jhe_vode*nhep*ehep
      rho_heat = uvb_density_a * (rho_vode/mean_rhob)**uvb_density_b
      heat = rho_heat*heat

      ! Convert back to code units
      ne_vode     = ne_vode / nh
      energy = (heat - lambda)*heat_from_cgs/(1.0d0+z_vode)**4

      ! Convert to the actual term to be used in e_out = e_in + dt*energy
      a = 1.d0 / (1.d0 + z_vode)
      energy = energy / rho_vode / a

end subroutine f_rhs

subroutine f_rhs_vec(time, e_in, energy)

      use amrex_fort_module, only : rt => amrex_real
      use fundamental_constants_module, only: e_to_cgs, density_to_cgs, & 
                                              heat_from_cgs
      use eos_module, only: iterate_ne_vec
      use atomic_rates_module, ONLY: tcoolmin, tcoolmax, ncooltab, deltat, &
                                     mproton, xhydrogen, &
                                     betah0, betahe0, betahep, betaff1, betaff4, &
                                     rechp, rechep, rechepp, &
                                     eh0, ehe0, ehep

      use amrex_error_module, only : amrex_abort
      use vode_aux_module       , only: t_vode_vec, ne_vode_vec, rho_vode_vec, z_vode
      use misc_params, only: simd_width

      implicit none

      real(rt),                        intent(in   ) :: time
      real(rt), dimension(simd_width), intent(inout) :: e_in
      real(rt), dimension(simd_width), intent(  out) :: energy

      real(rt), parameter :: compt_c = 1.01765467d-37, t_cmb = 2.725d0

      real(rt), dimension(simd_width) :: logT, tmp, fhi, flo
      real(rt), dimension(simd_width) :: ahp, ahep, ahepp, ad, geh0, gehe0, gehep
      real(rt), dimension(simd_width) :: bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp
      real(rt), dimension(simd_width) :: lambda_c, lambda_ff, lambda, heat
      real(rt), dimension(simd_width) :: rho, U
      real(rt) :: a
      real(rt), dimension(simd_width) :: nh, nh0, nhp, nhe0, nhep, nhepp
      integer, dimension(simd_width) :: j
      integer :: m
      logical, dimension(simd_width) :: hot

      do m = 1, simd_width
        if (e_in(m) .lt. 0.d0) then
           e_in(m) = tiny(e_in(m))
        endif
      end do

     ! Converts from code units to CGS
      rho = rho_vode_vec(1:simd_width) * density_to_cgs * (1.0d0+z_vode)**3
        u = e_in * e_to_cgs
      nh  = rho*xhydrogen/mproton

      if (time .gt. 1) then
         print *,'TIME INTO F_RHS ',time
         call amrex_abort("TOO BIG TIME IN F_RHS")
      end if

      ! Get gas temperature and individual ionization species
      call iterate_ne_vec(z_vode, u, t_vode_vec, nh, ne_vode_vec, nh0, nhp, nhe0, nhep, nhepp, simd_width)

      ! Convert species to CGS units: 
      ne_vode_vec(1:simd_width) = nh * ne_vode_vec(1:simd_width)
      nh0   = nh * nh0
      nhp   = nh * nhp
      nhe0  = nh * nhe0
      nhep  = nh * nhep
      nhepp = nh * nhepp

      logt = dlog10(t_vode_vec(1:simd_width))
      do m = 1, simd_width
         if (logt(m) .ge. tcoolmax) then ! Only free-free and Compton cooling are relevant
            lambda_ff(m) = 1.42d-27 * dsqrt(t_vode_vec(m)) * (1.1d0 + 0.34d0*dexp(-(5.5d0 - logt(m))**2 / 3.0d0)) &
                                 * (nhp(m) + 4.0d0*nhepp(m))*ne_vode_vec(m)
            lambda_c(m)  = compt_c*t_cmb**4 * ne_vode_vec(m) * (t_vode_vec(m) - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4

            energy(m)  = (-lambda_ff(m) -lambda_c(m)) * heat_from_cgs/(1.0d0+z_vode)**4

            ! Convert to the actual term to be used in e_out = e_in + dt*energy
            energy(m)  = energy(m) / rho_vode_vec(m) * (1.0d0+z_vode)
            ne_vode_vec(m) = ne_vode_vec(m) / nh(m)
            hot(m) = .true.
         else
            hot(m) = .false.
         endif
      end do

      do m = 1, simd_width
         if (.not. hot(m)) then
            ! Temperature floor
            if (logt(m) .le. tcoolmin) logt(m) = tcoolmin + 0.5d0*deltat
      
            ! Interpolate rates
            tmp(m) = (logt(m)-tcoolmin)/deltat
            j(m) = int(tmp(m))
            fhi(m) = tmp(m) - j(m)
            flo(m) = 1.0d0 - fhi(m)
            j(m) = j(m) + 1 ! F90 arrays start with 1
      
            bh0(m)   = flo(m)*betah0(j(m)) + fhi(m)*betah0(j(m)+1)
            bhe0(m)  = flo(m)*betahe0(j(m)) + fhi(m)*betahe0(j(m)+1)
            bhep(m)  = flo(m)*betahep(j(m)) + fhi(m)*betahep(j(m)+1)
            bff1(m)  = flo(m)*betaff1(j(m)) + fhi(m)*betaff1(j(m)+1)
            bff4(m)  = flo(m)*betaff4(j(m)) + fhi(m)*betaff4(j(m)+1)
            rhp(m)   = flo(m)*rechp(j(m)) + fhi(m)*rechp(j(m)+1)
            rhep(m)  = flo(m)*rechep(j(m)) + fhi(m)*rechep(j(m)+1)
            rhepp(m) = flo(m)*rechepp(j(m)) + fhi(m)*rechepp(j(m)+1)
      
            ! Cooling: 
            lambda(m) = ( bh0(m)*nh0(m) + bhe0(m)*nhe0(m) + bhep(m)*nhep(m) + &
                       rhp(m)*nhp(m) + rhep(m)*nhep(m) + rhepp(m)*nhepp(m) + &
                       bff1(m)*(nhp(m)+nhep(m)) + bff4(m)*nhepp(m) ) * ne_vode_vec(m)

            lambda_c(m) = compt_c*t_cmb**4*ne_vode_vec(m)*(t_vode_vec(m) - t_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4   ! Compton cooling
            lambda(m) = lambda(m) + lambda_c(m)
      
            ! Heating terms
            heat(m) = nh0(m)*eh0 + nhe0(m)*ehe0 + nhep(m)*ehep
      
            ! Convert back to code units
            ne_vode_vec(m)     = ne_vode_vec(m) / nh(m)
            energy(m) = (heat(m) - lambda(m))*heat_from_cgs/(1.0d0+z_vode)**4
      
            ! Convert to the actual term to be used in e_out = e_in + dt*energy
            a = 1.d0 / (1.d0 + z_vode)
            energy(m) = energy(m) / rho_vode_vec(m) / a
         end if
      end do

end subroutine f_rhs_vec


subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer         , intent(in   ) :: neq, ml, mu, nrpd, ipar
  real(rt), intent(in   ) :: y(neq), rpar, t
  real(rt), intent(  out) :: pd(neq,neq)

  ! Should never get here, we are using a numerical Jacobian
  print *,'IN JAC ROUTINE'
  stop

end subroutine jac
````

