
# File init\_e\_from\_rhoe\_3d.f90

[**File List**](files.md) **>** [**Initialization**](dir_71a4420ed1f8982e7234eb6a0b7e6d5d.md) **>** [**init\_e\_from\_rhoe\_3d.f90**](init__e__from__rhoe__3d_8f90.md)

[Go to the documentation of this file.](init__e__from__rhoe__3d_8f90.md) 


````cpp
      subroutine fort_init_e_from_rhoe(state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns, &
                                       lo,hi,a_old) bind(C, name="fort_init_e_from_rhoe")

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : urho, umx, umy, umz, ueint, ueden
      use  eos_params_module

      implicit none

      integer  ,intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns
      real(rt) ,intent(inout) :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      integer  ,intent(in   ) :: lo(3), hi(3)
      real(rt) ,intent(in   ) :: a_old

      ! Local variables
      integer          :: i,j,k
      !
      ! Compute energy from the EOS
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               state(i,j,k,ueden) = state(i,j,k,ueint) + &
                  0.5d0 * (state(i,j,k,umx)**2 + state(i,j,k,umy)**2 + state(i,j,k,umz)**2) / state(i,j,k,urho)

            enddo
         enddo
      enddo

      end subroutine fort_init_e_from_rhoe
````

