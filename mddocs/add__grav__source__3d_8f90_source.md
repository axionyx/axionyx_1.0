
# File add\_grav\_source\_3d.f90

[**File List**](files.md) **>** [**HydroFortran**](dir_1fab266cd447ad3f3624320661f845f1.md) **>** [**add\_grav\_source\_3d.f90**](add__grav__source__3d_8f90.md)

[Go to the documentation of this file.](add__grav__source__3d_8f90.md) 


````cpp
! :::
! ::: ------------------------------------------------------------------
! :::

    subroutine fort_add_grav_source(lo,hi,&
                               uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                               uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                               grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                               dt,a_old,a_new) &
                               bind(C, name="fort_add_grav_source")

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : nvar, urho, umx, umy, umz, &
           ueden, grav_source_type

      implicit none

      integer lo(3), hi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer  gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3

      real(rt)  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
      real(rt) uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      real(rt) grav(  gv_l1:  gv_h1,  gv_l2:  gv_h2,  gv_l3:  gv_h3,3)
      real(rt) dt
      real(rt) a_old, a_new

      real(rt) :: a_half, a_oldsq, a_newsq, a_newsq_inv
      real(rt) :: rho
      real(rt) :: SrU, SrV, SrW, SrE
      real(rt) :: rhoInv, dt_a_new
      real(rt) :: old_rhoeint, new_rhoeint, old_ke, new_ke
      integer          :: i, j, k

      a_half  = 0.5d0 * (a_old + a_new)
      a_oldsq = a_old * a_old
      a_newsq = a_new * a_new
      a_newsq_inv = 1.d0 / a_newsq

      dt_a_new    = dt / a_new

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 3: Puts all gravitational work into KE, not (rho e)

      ! Add gravitational source terms
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_ke = 0.5d0 * (uout(i,j,k,umx)**2 + uout(i,j,k,umy)**2 + uout(i,j,k,umz)**2) / &
                                 uout(i,j,k,urho) 
               old_rhoeint = uout(i,j,k,ueden) - old_ke
               ! ****   End Diagnostics ****

               rho    = uin(i,j,k,urho)
               rhoinv = 1.0d0 / rho

               sru = rho * grav(i,j,k,1)
               srv = rho * grav(i,j,k,2)
               srw = rho * grav(i,j,k,3)

               ! We use a_new here because we think of d/dt(a rho u) = ... + (rho g)
               uout(i,j,k,umx)   = uout(i,j,k,umx) + sru * dt_a_new
               uout(i,j,k,umy)   = uout(i,j,k,umy) + srv * dt_a_new
               uout(i,j,k,umz)   = uout(i,j,k,umz) + srw * dt_a_new

               if (grav_source_type .eq. 1) then

                   ! This does work (in 1-d)
                   ! Src = rho u dot g, evaluated with all quantities at t^n
                   sre = uin(i,j,k,umx) * grav(i,j,k,1) + &
                         uin(i,j,k,umy) * grav(i,j,k,2) + &
                         uin(i,j,k,umz) * grav(i,j,k,3)
                   uout(i,j,k,ueden) = (a_newsq*uout(i,j,k,ueden) + sre * (dt*a_half)) * a_newsq_inv

               else if (grav_source_type .eq. 3) then

                   new_ke = 0.5d0 * (uout(i,j,k,umx)**2 + uout(i,j,k,umy)**2 + uout(i,j,k,umz)**2) / &
                                     uout(i,j,k,urho) 
                   uout(i,j,k,ueden) = old_rhoeint + new_ke
               else 
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: bogus grav_source_type")
               end if

            enddo
         enddo
      enddo

      end subroutine fort_add_grav_source
````

