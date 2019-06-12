
# File integrate\_state\_fcvode\_3d.f90

[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**integrate\_state\_fcvode\_3d.f90**](integrate__state__fcvode__3d_8f90.md)

[Go to the documentation of this file.](integrate__state__fcvode__3d_8f90.md) 


````cpp
subroutine integrate_state_fcvode(lo, hi, &
                                  state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                  diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                  a, half_dt, min_iter, max_iter)
!
!   Calculates the sources to be added later on.
!
!   Parameters
!   ----------
!   lo : double array (3)
!       The low corner of the current box.
!   hi : double array (3)
!       The high corner of the current box.
!   state_* : double arrays
!       The state vars
!   diag_eos_* : double arrays
!       Temp and Ne
!   src_* : doubles arrays
!       The source terms to be added to state (iterative approx.)
!   double array (3)
!       The low corner of the entire domain
!   a : double
!       The current a
!   half_dt : double
!       time step size, in Mpc km^-1 s ~ 10^12 yr.
!
!   Returns
!   -------
!   state : double array (dims) @todo
!       The state vars
!
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : nvar, urho, ueden, ueint, &
                                   ndiag, temp_comp, ne_comp, zhi_comp, &
                                   gamma_minus_1
    use amrex_constants_module, only: m_pi
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_t_given_re, nyx_eos_given_rt
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_omb
    use comoving_nd_module, only: fort_integrate_comoving_a
    use atomic_rates_module, only: yhelium
    use vode_aux_module    , only: jh_vode, jhe_vode, z_vode, i_vode, j_vode, k_vode
    use reion_aux_module   , only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                   t_zhi, t_zheii, inhomogeneous_on
    use cvode_interface
    use fnvector_serial
    use fcvode_extras
    use, intrinsic :: iso_c_binding

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in)    :: a, half_dt
    integer         , intent(inout) :: max_iter, min_iter

    integer :: i, j, k
    real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
    real(rt) :: T_orig, ne_orig, e_orig
    real(rt) :: T_out , ne_out , e_out, mu, mean_rhob, T_H, T_He
    real(rt) :: species(5)

    integer(c_int) :: ierr       ! error flag from C functions
    real(c_double) :: tstart     ! initial time
    real(c_double) :: atol, rtol
    type(c_ptr) :: sunvec_y      ! sundials vector
    type(c_ptr) :: CVmem         ! CVODE memory
    integer(c_long), parameter :: neq = 1
    real(c_double), pointer :: yvec(:)

    allocate(yvec(neq))

    z = 1.d0/a - 1.d0
    call fort_integrate_comoving_a(a, a_end, half_dt)
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

    sunvec_y = n_vmake_serial(neq, yvec)
    if (.not. c_associated(sunvec_y)) then
        call amrex_abort('integrate_state_fcvode: sunvec = NULL')
    end if

    cvmem = fcvodecreate(cv_bdf, cv_newton)
    if (.not. c_associated(cvmem)) then
        call amrex_abort('integrate_state_fcvode: CVmem = NULL')
    end if

    tstart = 0.0
    ! CVodeMalloc allocates variables and initialize the solver. We can initialize the solver with junk because once we enter the
    ! (i,j,k) loop we will immediately call fcvreinit which reuses the same memory allocated from CVodeMalloc but sets up new
    ! initial conditions.
    ierr = fcvodeinit(cvmem, c_funloc(rhsfn), tstart, sunvec_y)
    if (ierr /= 0) then
       call amrex_abort('integrate_state_fcvode: FCVodeInit() failed')
    end if

    ! Set dummy tolerances. These will be overwritten as soon as we enter the loop and reinitialize the solver.
    rtol = 1.0d-5
    atol = 1.0d-10
    ierr = fcvodesstolerances(cvmem, rtol, atol)
    if (ierr /= 0) then
      call amrex_abort('integrate_state_fcvode: FCVodeSStolerances() failed')
    end if

    ierr = fcvdiag(cvmem)
    if (ierr /= 0) then
       call amrex_abort('integrate_state_fcvode: FCVDense() failed')
    end if

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)

                ! Original values
                rho     = state(i,j,k,urho)
                e_orig  = state(i,j,k,ueint) / rho
                t_orig  = diag_eos(i,j,k,temp_comp)
                ne_orig = diag_eos(i,j,k,  ne_comp)

                if (inhomogeneous_on) then
                   h_reion_z = diag_eos(i,j,k,zhi_comp)
                   if (z .gt. h_reion_z) then
                      jh_vode = 0
                   else
                      jh_vode = 1
                   endif
                endif

                if (e_orig .lt. 0.d0) then
                    !$OMP CRITICAL
                    print *,'negative e entering strang integration ',z, i,j,k, rho/mean_rhob, e_orig
                    call amrex_abort('bad e in strang')
                    !$OMP END CRITICAL
                end if

                i_vode = i
                j_vode = j
                k_vode = k

                call fcvode_wrapper(half_dt,rho,t_orig,ne_orig,e_orig,neq,cvmem,sunvec_y,yvec, &
                                              t_out ,ne_out ,e_out)

                if (e_out .lt. 0.d0) then
                    !$OMP CRITICAL
                    print *,'negative e exiting strang integration ',z, i,j,k, rho/mean_rhob, e_out
                    call flush(6)
                    !$OMP END CRITICAL
                    t_out  = 10.0
                    ne_out = 0.0
                    mu     = (1.0d0+4.0d0*yhelium) / (1.0d0+yhelium+ne_out)
                    e_out  = t_out / (gamma_minus_1 * mp_over_kb * mu)
!                    call amrex_abort('bad e out of strang')
                end if

                ! Update T and ne (do not use stuff computed in f_rhs, per vode manual)
                call nyx_eos_t_given_re(jh_vode, jhe_vode, t_out, ne_out, rho, e_out, a, species)

                ! Instanteneous heating from reionization:
                t_h = 0.0d0
                if (inhomogeneous_on .or. flash_h) then
                   if ((h_reion_z  .lt. z) .and. (h_reion_z  .ge. z_end)) t_h  = (1.0d0 - species(2))*max((t_zhi-t_out), 0.0d0)
                endif

                t_he = 0.0d0
                if (flash_he) then
                   if ((he_reion_z .lt. z) .and. (he_reion_z .ge. z_end)) t_he = (1.0d0 - species(5))*max((t_zheii-t_out), 0.0d0)
                endif

                if ((t_h .gt. 0.0d0) .or. (t_he .gt. 0.0d0)) then
                   t_out = t_out + t_h + t_he                            ! For simplicity, we assume
                   ne_out = 1.0d0 + yhelium                              !    completely ionized medium at
                   if (t_he .gt. 0.0d0) ne_out = ne_out + yhelium        !    this point.  It's a very minor
                   mu = (1.0d0+4.0d0*yhelium) / (1.0d0+yhelium+ne_out)   !    detail compared to the overall approximation.
                   e_out  = t_out / (gamma_minus_1 * mp_over_kb * mu)
                   call nyx_eos_t_given_re(jh_vode, jhe_vode, t_out, ne_out, rho, e_out, a, species)
                endif

                ! Update (rho e) and (rho E)
                state(i,j,k,ueint) = state(i,j,k,ueint) + rho * (e_out-e_orig)
                state(i,j,k,ueden) = state(i,j,k,ueden) + rho * (e_out-e_orig)

                ! Update T and ne
                diag_eos(i,j,k,temp_comp) = t_out
                diag_eos(i,j,k,  ne_comp) = ne_out

            end do ! i
        end do ! j
    end do ! k

    call n_vdestroy_serial(sunvec_y)
    call fcvodefree(cvmem)

    deallocate(yvec)

end subroutine integrate_state_fcvode
````

