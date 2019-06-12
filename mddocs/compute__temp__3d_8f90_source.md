
# File compute\_temp\_3d.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Src\_3d**](dir_723248e6e98dc7cb10ec13b7569a328c.md) **>** [**compute\_temp\_3d.f90**](compute__temp__3d_8f90.md)

[Go to the documentation of this file.](compute__temp__3d_8f90.md) 


````cpp

      ! *************************************************************************************

      subroutine fort_compute_temp(lo,hi, &
                                   state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                   diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                   comoving_a, print_fortran_warnings) &
      bind(C, name = "fort_compute_temp")

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use atomic_rates_module, only: this_z
      use meth_params_module, only : nvar, urho, umx, umy, umz, ueint, ueden, &
                                     ndiag, temp_comp, ne_comp, zhi_comp, &
                                     small_temp, heat_cool_type
      use reion_aux_module,    only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                     inhomogeneous_on
      use  eos_params_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer         , intent(in   ) :: print_fortran_warnings
      real(rt), intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(in   ) :: comoving_a

      integer          :: i,j,k, JH, JHe
      real(rt) :: rhoInv,eint
      real(rt) :: ke,dummy_pres
      real(rt) :: z

      z = 1.d0/comoving_a - 1.d0

      ! Flash reionization?
      if ((flash_h .eqv. .true.) .and. (z .gt. zhi_flash)) then
         jh = 0
      else
         jh = 1
      endif
      if ((flash_he .eqv. .true.) .and. (z .gt. zheii_flash)) then
         jhe = 0
      else
         jhe = 1
      endif

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (state(i,j,k,urho) <= 0.d0) then
                  print *,'   '
                  print *,'>>> Error: compute_temp ',i,j,k
                  print *,'>>> ... negative density ',state(i,j,k,urho)
                  print *,'    '
                  call amrex_error("Error:: compute_temp_3d.f90 :: compute_temp")
               end if
            enddo
         enddo
      enddo

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rhoinv = 1.d0 / state(i,j,k,urho)

               if (state(i,j,k,ueint) > 0.d0) then

                   eint = state(i,j,k,ueint) * rhoinv

                   jh = 1
                   if (inhomogeneous_on) then
                       if (z .gt. diag_eos(i,j,k,zhi_comp)) jh = 0
                   end if

                   call nyx_eos_t_given_re(jh, jhe, diag_eos(i,j,k,temp_comp), diag_eos(i,j,k,ne_comp), &
                                           state(i,j,k,urho), eint, comoving_a)

               else
                  if (print_fortran_warnings .gt. 0) then
                     print *,'   '
                     print *,'>>> Warning: (rho e) is negative in compute_temp: ',i,j,k
                  end if
                   ! Set temp to small_temp and compute corresponding internal energy
                   call nyx_eos_given_rt(eint, dummy_pres, state(i,j,k,urho), small_temp, &
                                         diag_eos(i,j,k,ne_comp), comoving_a)

                   ke = 0.5d0 * (state(i,j,k,umx)**2 + state(i,j,k,umy)**2 + state(i,j,k,umz)**2) * rhoinv

                   diag_eos(i,j,k,temp_comp) = small_temp
                   state(i,j,k,ueint) = state(i,j,k,urho) * eint
                   state(i,j,k,ueden) = state(i,j,k,ueint) + ke

               end if

            enddo
         enddo
      enddo

      end subroutine fort_compute_temp

      subroutine fort_compute_temp_vec(lo,hi, &
                                   state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                   diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                   comoving_a, print_fortran_warnings) &
      bind(C, name = "fort_compute_temp_vec")

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use atomic_rates_module, only: this_z
      use meth_params_module, only : nvar, urho, umx, umy, umz, ueint, ueden, &
                                     ndiag, temp_comp, ne_comp, small_temp, heat_cool_type
      use  eos_params_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer         , intent(in   ) :: print_fortran_warnings
      real(rt), intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(in   ) :: comoving_a

      integer          :: i,j,k
      real(rt) :: rhoInv,eint
      real(rt), dimension(hi(1)-lo(1)+1) :: ke,dummy_pres,small_temp_vec
      real(rt) :: z
      real(rt), dimension(hi(1)-lo(1)+1,4) :: eos_inputs_pos_ueint, eos_inputs_neg_ueint
      integer :: orig_indices(hi(1)-lo(1)+1,3)
      integer :: pos_eos_count, neg_eos_count

      z = 1.d0/comoving_a - 1.d0

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (state(i,j,k,urho) <= 0.d0) then
                  print *,'   '
                  print *,'>>> Error: compute_temp ',i,j,k
                  print *,'>>> ... negative density ',state(i,j,k,urho)
                  print *,'    '
                  call amrex_error("Error:: compute_temp_3d.f90 :: compute_temp")
               end if
            enddo
         enddo
      enddo

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)

            pos_eos_count = 0
            neg_eos_count = 0

            do i = lo(1),hi(1)
               rhoinv = 1.d0 / state(i,j,k,urho)

               if (state(i,j,k,ueint) > 0.d0) then

                   pos_eos_count = pos_eos_count + 1

                   eos_inputs_pos_ueint(pos_eos_count,1) = diag_eos(i,j,k,temp_comp)
                   eos_inputs_pos_ueint(pos_eos_count,2) = diag_eos(i,j,k,ne_comp)
                   eos_inputs_pos_ueint(pos_eos_count,3) = state(i,j,k,urho)
                   eos_inputs_pos_ueint(pos_eos_count,4) = state(i,j,k,ueint)*rhoinv

                   orig_indices(pos_eos_count,1) = i
                   orig_indices(pos_eos_count,2) = j
                   orig_indices(pos_eos_count,3) = k

               else

                   neg_eos_count = neg_eos_count + 1

                   eos_inputs_neg_ueint(neg_eos_count,1) = diag_eos(i,j,k,temp_comp) ! DON'T NEED THIS; GET RID OF IT
                   eos_inputs_neg_ueint(neg_eos_count,2) = diag_eos(i,j,k,ne_comp)
                   eos_inputs_neg_ueint(neg_eos_count,3) = state(i,j,k,urho)
                   eos_inputs_neg_ueint(neg_eos_count,4) = state(i,j,k,ueint)

                   orig_indices(neg_eos_count,1) = i
                   orig_indices(neg_eos_count,2) = j
                   orig_indices(neg_eos_count,3) = k

               end if
             end do

             ! For cells with positive E_int
             call nyx_eos_t_given_re_vec(eos_inputs_pos_ueint(1:pos_eos_count,1), &
                                         eos_inputs_pos_ueint(1:pos_eos_count,2), &
                                         eos_inputs_pos_ueint(1:pos_eos_count,3), &
                                         eos_inputs_pos_ueint(1:pos_eos_count,4), &
                                         comoving_a, &
                                         pos_eos_count)
             diag_eos(orig_indices(1:pos_eos_count,1),j,k,temp_comp) = eos_inputs_pos_ueint(1:pos_eos_count,1)
             diag_eos(orig_indices(1:pos_eos_count,1),j,k,ne_comp)   = eos_inputs_pos_ueint(1:pos_eos_count,2)

             ! For cells with negative E_int
             call nyx_eos_given_rt_vec(eos_inputs_neg_ueint(1:neg_eos_count,4), &
                                   dummy_pres(1:neg_eos_count), &
                                   eos_inputs_neg_ueint(1:neg_eos_count,3), &
                                   small_temp_vec(1:neg_eos_count), &
                                   eos_inputs_neg_ueint(1:neg_eos_count,2), &
                                   comoving_a, &
                                   neg_eos_count)

             ke(1:neg_eos_count) = 0.5d0 * (state(orig_indices(1:neg_eos_count,1),j,k,umx)*state(orig_indices(1:neg_eos_count,1),j,k,umx) + &
                                   state(orig_indices(1:neg_eos_count,1),j,k,umy)*state(orig_indices(1:neg_eos_count,1),j,k,umy) + &
                                   state(orig_indices(1:neg_eos_count,1),j,k,umz)*state(orig_indices(1:neg_eos_count,1),j,k,umz)) * rhoinv

             diag_eos(orig_indices(1:neg_eos_count,1),j,k,temp_comp) = small_temp_vec(1:neg_eos_count)
             state(orig_indices(1:neg_eos_count,1),j,k,ueint) = eos_inputs_neg_ueint(1:neg_eos_count,3) * eos_inputs_neg_ueint(1:neg_eos_count,4)
             state(orig_indices(1:neg_eos_count,1),j,k,ueden) = eos_inputs_neg_ueint(1:neg_eos_count,4) + ke(1:neg_eos_count)

         enddo
      enddo

      end subroutine fort_compute_temp_vec

      subroutine fort_compute_rho_temp(lo,hi,dx, &
                                     state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                                  diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                                  rho_ave,rho_T_sum, &
                                  T_sum,Tinv_sum,T_meanrho_sum,rho_sum,vol_sum,vol_mn_sum) &
      bind(C, name = "fort_compute_rho_temp")

      use meth_params_module, only : nvar, urho, ndiag, temp_comp

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt), intent(in   ) :: dx(3)
      real(rt), intent(in   ) :: rho_ave
      real(rt), intent(in   ) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(inout) :: rho_T_sum, rho_sum, T_sum, Tinv_sum, T_meanrho_sum
      real(rt), intent(inout) :: vol_sum, vol_mn_sum

      integer          :: i,j,k
      real(rt) :: rho_hi, rho_lo, vol

      vol = dx(1)*dx(2)*dx(3)
      rho_hi = 1.1d0*rho_ave
      rho_lo = 0.9d0*rho_ave
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                   t_sum =     t_sum + vol*diag_eos(i,j,k,temp_comp)
                tinv_sum =  tinv_sum + state(i,j,k,urho)/diag_eos(i,j,k,temp_comp)
               rho_t_sum = rho_t_sum + state(i,j,k,urho)*diag_eos(i,j,k,temp_comp)
                 rho_sum =   rho_sum + state(i,j,k,urho)
                 if ( (state(i,j,k,urho) .lt. rho_hi) .and. &
                      (state(i,j,k,urho) .gt. rho_lo) .and. &
                      (diag_eos(i,j,k,temp_comp) .le. 1.0e5) ) then
                         t_meanrho_sum = t_meanrho_sum + vol*dlog10(diag_eos(i,j,k,temp_comp))
                         vol_mn_sum = vol_mn_sum + vol
                 endif
                 vol_sum = vol_sum + vol
            enddo
         enddo
      enddo

      end subroutine fort_compute_rho_temp

      subroutine fort_compute_gas_frac(lo,hi,dx, &
                                       state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                                       diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                                       rho_ave, T_cut, rho_cut, &
                                       whim_mass, whim_vol, hh_mass, &
                                       hh_vol, igm_mass, igm_vol, mass_sum, vol_sum) &
      bind(C, name = "fort_compute_gas_frac")

      use meth_params_module, only : nvar, urho, ndiag, temp_comp

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt), intent(in   ) :: dx(3)
      real(rt), intent(in   ) :: rho_ave, T_cut, rho_cut
      real(rt), intent(in   ) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(inout) :: whim_mass, whim_vol, hh_mass, hh_vol, igm_mass, igm_vol
      real(rt), intent(inout) :: mass_sum, vol_sum

      integer :: i,j,k
      real(rt) :: vol, T, R

      vol = dx(1)*dx(2)*dx(3)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                 t = diag_eos(i,j,k,temp_comp)
                 r = state(i,j,k,urho) / rho_ave
                 if ( (t .gt. t_cut) .and. (r .le. rho_cut) ) then
                     whim_mass = whim_mass + state(i,j,k,urho)*vol
                     whim_vol  = whim_vol  + vol
                 else if ( (t .gt. t_cut) .and. (r .gt. rho_cut) ) then
                     hh_mass = hh_mass + state(i,j,k,urho)*vol
                     hh_vol  = hh_vol  + vol
                 else if ( (t .le. t_cut) .and. (r .le. rho_cut) ) then
                     igm_mass = igm_mass + state(i,j,k,urho)*vol
                     igm_vol  = igm_vol  + vol
                 endif
                 mass_sum = mass_sum + state(i,j,k,urho)*vol
                 vol_sum  = vol_sum + vol
            enddo
         enddo
      enddo

      end subroutine fort_compute_gas_frac

      subroutine fort_compute_max_temp_loc(lo,hi, &
                                           state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                           diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                           max_temp, den_maxt, imax, jmax, kmax) &
      bind(C, name = "fort_compute_max_temp_loc")

      use meth_params_module, only : temp_comp, nvar, urho, ndiag

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt), intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(in   ) :: max_temp
      real(rt), intent(  out) :: den_maxt
      integer         , intent(inout) :: imax,jmax,kmax

      integer                         :: i,j,k
      real(rt)                :: one_minus_eps

      one_minus_eps = 1.d0 - 1.d-12

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (diag_eos(i,j,k,temp_comp) .ge. one_minus_eps*max_temp) then
                  imax = i
                  jmax = j
                  kmax = k
                  den_maxt = state(i,j,k,urho)
               end if
            enddo
         enddo
      enddo

      end subroutine fort_compute_max_temp_loc
````

