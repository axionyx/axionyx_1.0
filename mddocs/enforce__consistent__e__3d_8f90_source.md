
# File enforce\_consistent\_e\_3d.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Src\_3d**](dir_723248e6e98dc7cb10ec13b7569a328c.md) **>** [**enforce\_consistent\_e\_3d.f90**](enforce__consistent__e__3d_8f90.md)

[Go to the documentation of this file.](enforce__consistent__e__3d_8f90.md) 


````cpp
   subroutine fort_enforce_consistent_e(lo,hi,state, &
                                        state_l1,state_l2,state_l3,state_h1,state_h2,state_h3) & 
     bind(C,name="fort_enforce_consistent_e")

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module, only : nvar, urho, umx, umy, umz, ueden, ueint

     implicit none

     integer          :: lo(3), hi(3)
     integer          :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     real(rt) :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     ! Local variables
     integer          :: i,j,k
     real(rt) :: u, v, w, rhoInv

     ! 
     ! Make sure to enforce (rho E) = (rho e) + 1/2 rho (u^2 +_ v^2 + w^2)
     !
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              rhoinv = 1.0d0 / state(i,j,k,urho)

              u = state(i,j,k,umx) * rhoinv
              v = state(i,j,k,umy) * rhoinv
              w = state(i,j,k,umz) * rhoinv

              state(i,j,k,ueden) = state(i,j,k,ueint) + &
                     0.5d0 * state(i,j,k,urho) * (u*u + v*v + w*w)

           end do
        end do
     end do

   end subroutine fort_enforce_consistent_e
````

