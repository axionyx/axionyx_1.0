
# File EstDt\_3d.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Src\_3d**](dir_723248e6e98dc7cb10ec13b7569a328c.md) **>** [**EstDt\_3d.f90**](EstDt__3d_8f90.md)

[Go to the documentation of this file.](EstDt__3d_8f90.md) 


````cpp
     subroutine fort_estdt(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,dx,dt,a_old) &
        bind(C, name = "fort_estdt")

     use amrex_fort_module, only : rt => amrex_real
     use eos_module
     use meth_params_module, only : nvar, urho, umx, umy, umz, ueint
     use  eos_params_module

     implicit none
     ! 
     ! NOTE: for comoving coordinates, the factor of "a" is multiplied *outside* this routine
     ! 
     integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
     integer          :: lo(3), hi(3)
     real(rt) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
     real(rt) :: dx(3), dt
     real(rt) :: a_old

     real(rt) :: e, c
     real(rt) :: rhoInv,ux,uy,uz,dt1,dt2,dt3
     real(rt) :: sqrtK,grid_scl,dt4
     integer          :: i,j,k

     real(rt), parameter :: onethird = 1.d0/3.d0

     grid_scl = (dx(1)*dx(2)*dx(3))**onethird
     !
     ! Translate to primitive variables, compute sound speed
     !
     do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rhoinv = 1.d0 / u(i,j,k,urho)
               ux     = u(i,j,k,umx)*rhoinv
               uy     = u(i,j,k,umy)*rhoinv
               uz     = u(i,j,k,umz)*rhoinv

               ! Use internal energy for calculating dt 
               e  = u(i,j,k,ueint)*rhoinv

               ! Protect against negative e
               if (e .gt. 0.d0) then
                  call nyx_eos_soundspeed(c,u(i,j,k,urho),e)
               else
                  c = 0.d0
               end if

               dt1 = dx(1)/(c + abs(ux))
               dt2 = dx(2)/(c + abs(uy))
               dt3 = dx(3)/(c + abs(uz))
               dt  = min(dt,dt1,dt2,dt3)
            enddo
         enddo
     enddo

     end subroutine fort_estdt
````

