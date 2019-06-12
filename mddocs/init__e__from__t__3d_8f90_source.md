
# File init\_e\_from\_t\_3d.f90

[**File List**](files.md) **>** [**Initialization**](dir_71a4420ed1f8982e7234eb6a0b7e6d5d.md) **>** [**init\_e\_from\_t\_3d.f90**](init__e__from__t__3d_8f90.md)

[Go to the documentation of this file.](init__e__from__t__3d_8f90.md) 


````cpp
      subroutine fort_init_e_from_t(state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns, &
                                    diag,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3,nd, &
                                   lo,hi,a_old) bind(C, name="fort_init_e_from_t")


      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : urho, umx, umy, umz, ueint, ueden, temp_comp, ne_comp
      use  eos_params_module

      implicit none

      integer , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns
      integer , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3,nd
      real(rt), intent(inout) :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      real(rt), intent(in   ) ::  diag(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,nd)
      integer , intent(in   ) :: lo(3), hi(3)
      real(rt), intent(in   ) :: a_old

      ! Local variables
      real(rt) :: e, pres, T
      integer  :: i,j,k
      !
      ! Compute energy from the EOS
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               t  = diag(i,j,k,temp_comp)

               ! Call EOS to get the internal energy
               call nyx_eos_given_rt(e, pres, state(i,j,k,urho), diag(i,j,k,temp_comp), &
                                     diag(i,j,k,ne_comp), a_old)

               state(i,j,k,ueint) = state(i,j,k,urho) * e

               state(i,j,k,ueden) = state(i,j,k,ueint) + &
                  0.5d0 * (state(i,j,k,umx)**2 + state(i,j,k,umy)**2 + state(i,j,k,umz)**2) / state(i,j,k,urho)

            enddo
         enddo
      enddo

      end subroutine fort_init_e_from_t
````

