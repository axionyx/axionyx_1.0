
# File atomic\_rates.f90

[**File List**](files.md) **>** [**EOS**](dir_2a6406f09975eea078703cc63b0e3416.md) **>** [**atomic\_rates.f90**](atomic__rates_8f90.md)

[Go to the documentation of this file.](atomic__rates_8f90.md) 


````cpp
! *************************************************************************************
! Tabulates cooling and UV heating rates. 
!
!     Units are CGS, temperature is in K
!
!     Two sets of rates, as in:
!      1) Katz, Weinberg & Hernquist, 1996: Astrophysical Journal Supplement v. 105, p.19
!      2) Lukic et al. 2015: Monthly Notices of the Royal Astronomical Society, v. 446, p.3697
! 
!     NOTE: This is executed only once per run, and rates are ugly, thus 
!           execution efficiency is not important, but readability of 
!           the code is. -- Zarija
!
! *************************************************************************************

module atomic_rates_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  ! Photo- rates (from file)
  integer, private :: NCOOLFILE
  real(rt), dimension(:), allocatable, private :: lzr
  real(rt), dimension(:), allocatable, private :: rggh0, rgghe0, rgghep
  real(rt), dimension(:), allocatable, private :: reh0, rehe0, rehep

  ! Other rates (from equations)
  integer, parameter, public :: NCOOLTAB=2000
  real(rt), dimension(NCOOLTAB+1), public :: AlphaHp, AlphaHep, AlphaHepp, Alphad
  real(rt), dimension(NCOOLTAB+1), public :: GammaeH0, GammaeHe0, GammaeHep
  real(rt), dimension(NCOOLTAB+1), public :: BetaH0, BetaHe0, BetaHep, Betaff1, Betaff4
  real(rt), dimension(NCOOLTAB+1), public :: RecHp, RecHep, RecHepp

  real(rt), public, save :: this_z, ggh0, gghe0, gghep, eh0, ehe0, ehep
 
  real(rt), parameter, public :: TCOOLMIN = 0.0d0, tcoolmax = 9.0d0  ! in log10
  real(rt), parameter, public :: TCOOLMIN_R = 10.0d0**tcoolmin, tcoolmax_r = 10.0d0**tcoolmax
  real(rt), parameter, public :: deltaT = (tcoolmax - tcoolmin)/ncooltab

  real(rt), parameter, public :: MPROTON = 1.6726231d-24, boltzmann = 1.3806e-16

  real(rt), public, save :: uvb_density_A = 1.0d0, uvb_density_b = 0.0d0, mean_rhob

  ! Note that XHYDROGEN can be set by a call to set_xhydrogen which now
  ! lives in set_method_params.
  real(rt), public :: XHYDROGEN = 0.76d0
  real(rt), public :: YHELIUM   = 7.8947368421d-2  ! (1.0d0-XHYDROGEN)/(4.0d0*XHYDROGEN)

  contains

      subroutine fort_tabulate_rates() bind(C, name='fort_tabulate_rates')
      use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
      use amrex_parmparse_module
      use amrex_constants_module, only: m_pi
      use fundamental_constants_module, only: gconst
      use comoving_module, only: comoving_h,comoving_omb
      use reion_aux_module, only: zhi_flash, zheii_flash, t_zhi, t_zheii, &
                                  flash_h, flash_he, inhomogeneous_on

      integer :: i, inhomo_reion
      logical, parameter :: Katz96=.false.
      real(rt), parameter :: t3=1.0d3, t5=1.0d5, t6=1.0d6
      real(rt) :: t, U, E, y, sqrt_t, corr_term, tmp
      logical, save :: first=.true.

      character(len=:), allocatable :: file_in
      type(amrex_parmparse) :: pp

      if (first) then

         first = .false.

         ! Get info from inputs
         call amrex_parmparse_build(pp, "nyx")
         call pp%query("inhomo_reion"             , inhomo_reion)
         call pp%query("uvb_rates_file"           , file_in)
         call pp%query("uvb_density_A"            , uvb_density_a)
         call pp%query("uvb_density_B"            , uvb_density_b)
         call pp%query("reionization_zHI_flash"   , zhi_flash)
         call pp%query("reionization_zHeII_flash" , zheii_flash)
         call pp%query("reionization_T_zHI"       , t_zhi)
         call pp%query("reionization_T_zHeII"     , t_zheii)
         call amrex_parmparse_destroy(pp)

         if (amrex_pd_ioprocessor()) then
            print*, 'TABULATE_RATES: reionization parameters are:'
            print*, '    reionization_zHI_flash     = ', zhi_flash
            print*, '    reionization_zHeII_flash   = ', zheii_flash
            print*, '    reionization_T_zHI         = ', t_zhi
            print*, '    reionization_T_zHeII       = ', t_zheii

            print*, 'TABULATE_RATES: rho-dependent heating parameters are:'
            print*, '    A       = ', uvb_density_a
            print*, '    B       = ', uvb_density_b
            print*, '    UVB heating rates will be multiplied by A*(rho/rho_mean)**B'
        endif

        ! Save mean density (in code units) for density-dependent heating
        mean_rhob = comoving_omb * 3.d0*(comoving_h*100.d0)**2 / (8.d0*m_pi*gconst)

         ! Set options in reion_aux_module
         !   Hydrogen reionization
         if (zhi_flash .gt. 0.0) then
            if (inhomo_reion .gt. 0) then
               if (amrex_pd_ioprocessor()) print*, 'TABULATE_RATES: ignoring reionization_zHI, as nyx.inhomo_reion > 0'
               flash_h = .false.
               inhomogeneous_on = .true.
            else
               flash_h = .true.
               inhomogeneous_on = .false.
            endif
         else
            flash_h = .false.
            if (inhomo_reion .gt. 0) then
               inhomogeneous_on = .true.
            else
               inhomogeneous_on = .false.
            endif
         endif

         !   Helium reionization
         if (zheii_flash .gt. 0.0) then
            flash_he = .true.
         else
            flash_he = .false.
         endif

         if (amrex_pd_ioprocessor()) then
            print*, 'TABULATE_RATES: reionization flags are set to:'
            print*, '    Hydrogen flash            = ', flash_h
            print*, '    Helium   flash            = ', flash_he
            print*, '    inhomogeneous_on (H only) = ', inhomogeneous_on
         endif


         ! Read in UVB rates from a file
         if (len(file_in) .gt. 0) then
            open(unit=11, file=file_in, status='old')
            if (amrex_pd_ioprocessor()) then
               print*, 'TABULATE_RATES: UVB file is set in inputs ('//file_in//').'
            endif
         else
            open(unit=11, file='TREECOOL', status='old')
            if (amrex_pd_ioprocessor()) then
               print*, 'TABULATE_RATES: UVB file is defaulted to "TREECOOL".'
            endif
         endif

         ncoolfile = 0
         do
            read(11,*,end=10) tmp, tmp, tmp, tmp, tmp,  tmp, tmp
            ncoolfile = ncoolfile + 1
         end do
         10 rewind(11)

         allocate( lzr(ncoolfile), rggh0(ncoolfile), rgghe0(ncoolfile), rgghep(ncoolfile) )
         allocate( reh0(ncoolfile), rehe0(ncoolfile), rehep(ncoolfile) )

         do i = 1, ncoolfile
            read(11,*) lzr(i), rggh0(i), rgghe0(i), rgghep(i), &
                                reh0(i),  rehe0(i),  rehep(i)
         end do
         close(11)

         ! Initialize cooling tables
         t = 10.0d0**tcoolmin
         if (katz96) then
            do i = 1, ncooltab+1
               ! Rates are as in Katz et al. 1996
               sqrt_t = dsqrt(t)      
               corr_term    = 1.d0 / (1.0d0 + sqrt_t/dsqrt(t5))

               ! Recombination rates
               ! Alphad: dielectronic recombination rate of singly ioniozed helium
               alphad(i)    = 1.90d-03/(t*sqrt_t) * dexp(-4.7d5/t) * (1.0d0+0.3d0*dexp(-9.4d4/t))
               alphahp(i)   = 8.40d-11/sqrt_t * (t/t3)**(-0.2d0) / (1.0d0 + (t/t6)**0.7d0)
               alphahep(i)  = 1.50d-10 * t**(-0.6353d0)
               alphahepp(i) = 3.36d-10/sqrt_t * (t/t3)**(-0.2d0) / (1.0d0 + (t/t6)**0.7d0)

               ! Collisional ionization rates
               gammaeh0(i)  = 5.85d-11*sqrt_t * dexp(-157809.1d0/t) * corr_term
               gammaehe0(i) = 2.38d-11*sqrt_t * dexp(-285335.4d0/t) * corr_term
               gammaehep(i) = 5.68d-12*sqrt_t * dexp(-631515.0d0/t) * corr_term

               ! Collisional ionization & excitation cooling rates
               betah0(i)  = 7.5d-19 * dexp(-118348.0d0/t) * corr_term + 2.171d-11*gammaeh0(i)
               betahe0(i) = 3.941d-11 * gammaehe0(i)
               betahep(i) = 5.54d-17 * t**(-0.397d0) * dexp(-473638.0d0/t) * corr_term + &
                            8.715d-11 * gammaehep(i)

               ! Recombination cooling rates
               rechp(i)   = 1.036d-16 * t * alphahp(i)
               rechep(i)  = 1.036d-16 * t * alphahep(i) + 6.526d-11 * alphad(i)
               rechepp(i) = 1.036d-16 * t * alphahepp(i)

               ! Free-free cooling rate
               betaff1(i) = 1.42d-27 * sqrt_t * (1.1d0 + 0.34d0*dexp(-(5.5d0 - dlog10(t))**2 / 3.0d0))
               betaff4(i) = betaff1(i)

               t = t*10.0d0**deltat
            enddo
         else
            do i = 1, ncooltab+1
               ! Rates are as in Lukic et al.
               sqrt_t = dsqrt(t)

               ! Recombination rates
               ! Alphad: dielectronic recombination rate of singly ioniozed helium
               alphad(i)    = 1.90d-03/(t*sqrt_t) * dexp(-4.7d5/t) * (1.0d0+0.3d0*dexp(-9.4d4/t))
               alphahp(i)   = 7.982d-11 / (dsqrt(t/3.148d0)* (1.0d0+dsqrt(t/3.148d0))**0.252 * &
                                           (1.0d0+dsqrt(t/7.036d5))**1.748)
               if (t .le. 1.0d6) then
                  alphahep(i)  = 3.294d-11 / (dsqrt(t/15.54d0)* (1.0d0+dsqrt(t/15.54d0))**0.309 * &
                                              (1.0d0+dsqrt(t/3.676d7))**1.691)
               else
                  alphahep(i)  = 9.356d-10 / (dsqrt(t/4.266d-2)* (1.0d0+dsqrt(t/4.266d-2))**0.2108 * &
                                              (1.0d0+dsqrt(t/4.677d6))**1.7892)
               endif
               alphahepp(i) = 1.891d-10 / (dsqrt(t/9.37d0)* (1.0d0+dsqrt(t/9.37d0))**0.2476 * &
                                           (1.0d0+dsqrt(t/2.774d6))**1.7524)

               ! Collisional ionization rates
               e = 13.6d0
               u = 1.16045d4*e/t
               gammaeh0(i)  = 2.91d-8*u**0.39*dexp(-u) / (0.232d0+u)
               e = 24.6d0
               u = 1.16045d4*e/t
               gammaehe0(i) = 1.75d-8*u**0.35*dexp(-u) / (0.18d0+u)
               e = 54.4d0
               u = 1.16045d4*e/t
               gammaehep(i) = 2.05d-9*(1.0d0+dsqrt(u))*u**0.25*dexp(-u) / (0.265d0+u)

               ! Collisional ionization & excitation cooling rates
               corr_term = 1.d0 / (1.0d0 + sqrt_t/dsqrt(5.0d7))
               y = dlog(t)
               if (t .le. 1.0d5) then
                  betah0(i)  = 1.0d-20 * dexp( 2.137913d2 - 1.139492d2*y + 2.506062d1*y**2 - &
                                               2.762755d0*y**3 + 1.515352d-1*y**4 - &
                                               3.290382d-3*y**5 - 1.18415d5/t )
               else
                  betah0(i)  = 1.0d-20 * dexp( 2.7125446d2 - 9.8019455d1*y + 1.400728d1*y**2 - &
                                               9.780842d-1*y**3 + 3.356289d-2*y**4 - &
                                               4.553323d-4*y**5 - 1.18415d5/t )
               endif
               betahe0(i) = 9.38d-22 * sqrt_t * dexp(-285335.4d0/t) * corr_term
               betahep(i) = (5.54d-17 * t**(-0.397d0) * dexp(-473638.0d0/t) + & 
                             4.95d-22 * sqrt_t * dexp(-631515.0d0/t) )*corr_term

               ! Recombination cooling rates
               rechp(i)   = 2.851d-27 * sqrt_t * (5.914d0-0.5d0*dlog(t)+1.184d-2*t**(1.0d0/3.0d0))
               rechep(i)  = 1.55d-26 * t**0.3647 + 1.24d-13/(t*sqrt_t) * dexp(-4.7d5/t) * & 
                                                                         (1.0d0+0.3d0*dexp(-9.4d4/t))
               rechepp(i) = 1.14d-26 * sqrt_t * (6.607d0-0.5d0*dlog(t)+7.459d-3*t**(1.0d0/3.0d0))

               ! Free-free cooling rate
               if (t .le. 3.2d5) then
                  betaff1(i) = 1.426d-27 * sqrt_t * (0.79464d0 + 0.1243d0*dlog10(t))
               else
                  betaff1(i) = 1.426d-27 * sqrt_t * (2.13164d0 - 0.1240d0*dlog10(t))
               endif

               if (t/4.0d0 .le. 3.2d5) then
                  betaff4(i) = 1.426d-27 * sqrt_t * 4.0d0*(0.79464d0 + 0.1243d0*dlog10(t/4.0d0))
               else
                  betaff4(i) = 1.426d-27 * sqrt_t * 4.0d0*(2.13164d0 - 0.1240d0*dlog10(t/4.0d0))
               endif
               
               t = t*10.0d0**deltat
            enddo
         endif  ! Katz rates

      end if  ! first_call

      end subroutine fort_tabulate_rates

      ! ****************************************************************************

      subroutine fort_interp_to_this_z(z) bind(C, name='fort_interp_to_this_z')

      use vode_aux_module, only: z_vode

      real(rt), intent(in) :: z
      real(rt) :: lopz, fact
      integer :: i, j

      this_z = z
      z_vode = z
      lopz   = dlog10(1.0d0 + z)

      if (lopz .ge. lzr(ncoolfile)) then
         ggh0  = 0.0d0
         gghe0 = 0.0d0
         gghep = 0.0d0
         eh0   = 0.0d0
         ehe0  = 0.0d0
         ehep  = 0.0d0
         return
      endif

      if (lopz .le. lzr(1)) then
         j = 1
      else
         do i = 2, ncoolfile
            if (lopz .lt. lzr(i)) then
               j = i-1
               exit
            endif
         enddo
      endif

      fact  = (lopz-lzr(j))/(lzr(j+1)-lzr(j))

      ggh0  = rggh0(j)  + (rggh0(j+1)-rggh0(j))*fact
      gghe0 = rgghe0(j) + (rgghe0(j+1)-rgghe0(j))*fact
      gghep = rgghep(j) + (rgghep(j+1)-rgghep(j))*fact
      eh0   = reh0(j)   + (reh0(j+1)-reh0(j))*fact
      ehe0  = rehe0(j)  + (rehe0(j+1)-rehe0(j))*fact
      ehep  = rehep(j)  + (rehep(j+1)-rehep(j))*fact

      end subroutine fort_interp_to_this_z

end module atomic_rates_module
````

