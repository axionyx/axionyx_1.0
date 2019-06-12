
# File Nyx\_grav\_sources\_3d.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**SourceTerms**](dir_7c1c0d2e2a0285e12a54f57a60f809aa.md) **>** [**Nyx\_grav\_sources\_3d.f90**](Nyx__grav__sources__3d_8f90.md)

[Go to the documentation of this file.](Nyx__grav__sources__3d_8f90.md) 


````cpp

! :::
! ::: ------------------------------------------------------------------
! :::

      !===========================================================================
      ! This is called from the C++ so the threading happens here...
      !===========================================================================
      subroutine fort_correct_gsrc(lo,hi, &
                              gold,gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3, &
                              gnew,gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3, &
                              uold,uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3, &
                              unew,unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3, &
                              a_old,a_new,dt) &
                              bind(C, name="fort_correct_gsrc")

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : nvar, urho, umx, umy, umz, ueden, grav_source_type

      implicit none

      integer lo(3),hi(3)
      integer gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3
      integer gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3
      integer uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3
      integer unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3
      real(rt)   gold(gold_l1:gold_h1,gold_l2:gold_h2,gold_l3:gold_h3,3)
      real(rt)   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,gnew_l3:gnew_h3,3)
      real(rt)  uold(uold_l1:uold_h1,uold_l2:uold_h2,uold_l3:uold_h3,NVAR)
      real(rt)  unew(unew_l1:unew_h1,unew_l2:unew_h2,unew_l3:unew_h3,NVAR)
      real(rt)  a_old,a_new,dt

      integer i,j,k
      real(rt) SrU_old, SrV_old, SrW_old
      real(rt) SrU_new, SrV_new, SrW_new
      real(rt) SrUcorr, SrVcorr, SrWcorr, SrEcorr
      real(rt) rhoo, Upo, Vpo, Wpo
      real(rt) rhon, Upn, Vpn, Wpn

      real(rt) a_half, a_newsq, rhooinv, rhoninv, a_new_inv
      real(rt) old_ke, old_rhoeint
      real(rt) new_ke, new_rhoeint

      a_half    = 0.5d0 * (a_old + a_new)
      a_newsq   = a_new*a_new
      a_new_inv = 1.0d0 / a_new

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 3: Puts all gravitational work into KE, not (rho e)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_ke = 0.5d0 * (unew(i,j,k,umx)**2 + unew(i,j,k,umy)**2 + unew(i,j,k,umz)**2) / &
                                 unew(i,j,k,urho) 
               old_rhoeint = unew(i,j,k,ueden) - old_ke
               ! ****   End Diagnostics ****

               rhoo    = uold(i,j,k,urho)
               rhooinv = 1.0d0 / uold(i,j,k,urho)
               upo     = uold(i,j,k,umx) * rhooinv
               vpo     = uold(i,j,k,umy) * rhooinv
               wpo     = uold(i,j,k,umz) * rhooinv

               ! Define old source terms
               sru_old = rhoo * gold(i,j,k,1)
               srv_old = rhoo * gold(i,j,k,2)
               srw_old = rhoo * gold(i,j,k,3)

               rhon    = unew(i,j,k,urho)
               rhoninv = 1.0d0 / unew(i,j,k,urho)
               upn     = unew(i,j,k,umx) * rhoninv
               vpn     = unew(i,j,k,umy) * rhoninv
               wpn     = unew(i,j,k,umz) * rhoninv

               ! Define new source terms
               sru_new = rhon * gnew(i,j,k,1)
               srv_new = rhon * gnew(i,j,k,2)
               srw_new = rhon * gnew(i,j,k,3)

               ! Define corrections to source terms
               srucorr = 0.5d0*(sru_new - sru_old)
               srvcorr = 0.5d0*(srv_new - srv_old)
               srwcorr = 0.5d0*(srw_new - srw_old)

               ! This does work (in 1-d)
               if (grav_source_type .eq. 1) then
                   srecorr =  0.5d0 * ( (sru_new * upn - sru_old * upo) + &
                                        (srv_new * vpn - srv_old * vpo) + &
                                        (srw_new * wpn - srw_old * wpo) )
               end if

               ! Correct state with correction terms
               unew(i,j,k,umx)   = unew(i,j,k,umx)   + srucorr*dt * a_new_inv
               unew(i,j,k,umy)   = unew(i,j,k,umy)   + srvcorr*dt * a_new_inv
               unew(i,j,k,umz)   = unew(i,j,k,umz)   + srwcorr*dt * a_new_inv

               if (grav_source_type .eq. 1) then
                   unew(i,j,k,ueden) = unew(i,j,k,ueden) + srecorr*dt * (a_half / a_newsq)
               else if (grav_source_type .eq. 3) then
                   new_ke = 0.5d0 * (unew(i,j,k,umx)**2 + unew(i,j,k,umy)**2 + unew(i,j,k,umz)**2) / &
                                     unew(i,j,k,urho) 
                   unew(i,j,k,ueden) = old_rhoeint + new_ke
               else 
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: bogus grav_source_type")
               end if

            enddo
         enddo
      enddo

      end subroutine fort_correct_gsrc
! :::
! ::: ------------------------------------------------------------------
! :::
      subroutine fort_syncgsrc(lo,hi, &
                              gphi,gphi_l1,gphi_l2,gphi_l3,gphi_h1,gphi_h2,gphi_h3, &
                              gdphi,gdphi_l1,gdphi_l2,gdphi_l3,gdphi_h1,gdphi_h2,gdphi_h3, &
                              state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                              dstate,dstate_l1,dstate_l2,dstate_l3, &
                              dstate_h1,dstate_h2,dstate_h3, &
                              sync_src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,a_new,dt) &
                              bind(C, name="fort_syncgsrc")

      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : nvar, urho, umx, umy, umz

      implicit none
 
      integer lo(3),hi(3)
      integer gphi_l1,gphi_l2,gphi_l3,gphi_h1,gphi_h2,gphi_h3
      integer gdphi_l1,gdphi_l2,gdphi_l3,gdphi_h1,gdphi_h2,gdphi_h3
      integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      integer dstate_l1,dstate_l2,dstate_l3,dstate_h1,dstate_h2,dstate_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      real(rt)   gphi(gphi_l1:gphi_h1,gphi_l2:gphi_h2,gphi_l3:gphi_h3,3)
      real(rt)  gdphi(gdphi_l1:gdphi_h1,gdphi_l2:gdphi_h2,gdphi_l3:gdphi_h3,3)
      real(rt)  state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
      real(rt) dstate(dstate_l1:dstate_h1,dstate_l2:dstate_h2,dstate_l3:dstate_h3,3+1)
      real(rt) sync_src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,3+1)
      real(rt) a_new,dt
 
      !    Note that dstate is drho and drhoU, state is the entire state, and src
      !    is S_rhoU and S_rhoE
 
      integer          :: i,j,k
      real(rt) :: rho_pre, rhoU_pre, rhoV_pre, rhoW_pre
      real(rt) :: gx, gy, gz, dgx, dgy, dgz, SrU, SrV, SrW, SrE, a_new_inv
 
      a_new_inv = 1.0d0 / a_new

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rho_pre  = state(i,j,k,urho) - dstate(i,j,k,1)
               rhou_pre = state(i,j,k,umx)  - dstate(i,j,k,2)
               rhov_pre = state(i,j,k,umy)  - dstate(i,j,k,3)
               rhow_pre = state(i,j,k,umz)  - dstate(i,j,k,4)

               gx  = gphi(i,j,k,1)
               gy  = gphi(i,j,k,2)
               gz  = gphi(i,j,k,3)

               dgx = gdphi(i,j,k,1)
               dgy = gdphi(i,j,k,2)
               dgz = gdphi(i,j,k,3)

               sru = dstate(i,j,k,1)*gx + rho_pre*dgx
               srv = dstate(i,j,k,1)*gy + rho_pre*dgy
               srw = dstate(i,j,k,1)*gz + rho_pre*dgz

               sre = ( sru * (rhou_pre + (0.5d0*dt)*sru) + &
                       srv * (rhov_pre + (0.5d0*dt)*srv) + &
                       srw * (rhow_pre + (0.5d0*dt)*srw) ) / rho_pre

               sync_src(i,j,k,1) = sru * a_new_inv
               sync_src(i,j,k,2) = srv * a_new_inv
               sync_src(i,j,k,3) = srw * a_new_inv
               sync_src(i,j,k,4) = sre * a_new_inv

            enddo
         enddo
      enddo

      end subroutine fort_syncgsrc
````

