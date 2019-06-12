
# File fcvode\_extras\_src.f90

[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**fcvode\_extras\_src.f90**](fcvode__extras__src_8f90.md)

[Go to the documentation of this file.](fcvode__extras__src_8f90.md) 


````cpp
module fcvode_extras_src

  implicit none

  contains

subroutine fcvode_wrapper_with_source(dt, rho_in, T_in, ne_in, e_in, neq, cvmem, &
                              sunvec_y, yvec, rho_out, T_out, ne_out, e_out, rho_src, e_src)

    use amrex_fort_module, only : rt => amrex_real
    use vode_aux_module, only: rho_vode, t_vode, ne_vode, &
                               i_vode, j_vode, k_vode, nr_vode, rho_src_vode, e_src_vode,&
                               rho_init_vode

    use eos_params_module
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_omb

    use cvode_interface
    use fnvector_serial
    use eos_module, only: vode_rtol, vode_atol_scaled
    use, intrinsic :: iso_c_binding

    implicit none

    real(rt), intent(in   ) :: dt
    real(rt), intent(in   ) :: rho_in, T_in, ne_in, e_in
    type(c_ptr), value :: cvmem
    type(c_ptr), value :: sunvec_y
    real(rt), intent(  out) ::  rho_out, T_out,ne_out,e_out
    real(rt),  intent(in   ) ::         rho_src, e_src

    real(rt) mean_rhob

    real(c_double), pointer, dimension(:) :: atol
    real(c_double) :: rtol
    real(c_double) :: time, tout
    integer(c_long), intent(in) :: neq
    real(c_double), pointer, intent(in) :: yvec(:)
    
    integer(c_int) :: ierr
    real(c_double) :: t_soln
    type(c_ptr) :: sunvec_atol
    logical, save :: firstCall = .true.

    allocate(atol(neq))

    t_vode   = t_in
    ne_vode  = ne_in
    rho_vode = rho_in
    rho_init_vode = rho_in
    nr_vode  = 0
    rho_src_vode = rho_src
    e_src_vode = e_src

    ! Initialize the integration time
    time = 0.d0
    
    ! We will integrate "e" and "rho" in time. 
    yvec(1) = e_in
    yvec(2) = rho_in

    ! Set the tolerances.  
    atol(1) = 1.d-4 * e_in
    atol(2) = 1.d-4 * rho_in
    rtol = 1.d-4

    sunvec_atol = n_vmake_serial(neq, atol)
    ierr = fcvodereinit(cvmem, time, sunvec_y)
    ierr = fcvodesvtolerances(cvmem, rtol, sunvec_atol)
    
    ierr = fcvode(cvmem, dt, sunvec_y, time, cv_normal)
    
    e_out  = yvec(1)
    rho_out = yvec(2)
!    rho_out = rho_init_vode + dt * rho_src_vode
    t_out  = t_vode
    ne_out = ne_vode

    if (ierr .ne. 0) then
       print *, 'istate = ', ierr, 'at (i,j,k) ',i_vode,j_vode,k_vode
       call amrex_error("ERROR in fcvode_wrapper_with_src: integration failed")
    endif

    call n_vdestroy_serial(sunvec_atol)

end subroutine fcvode_wrapper_with_source

subroutine fcvode_wrapper_with_source_single(dt, rho_in, T_in, ne_in, e_in, neq, cvmem, &
                              sunvec_y, yvec, rho_out, T_out, ne_out, e_out, rho_src, e_src)

    use amrex_fort_module, only : rt => amrex_real
    use vode_aux_module, only: rho_vode, t_vode, ne_vode, &
                               i_vode, j_vode, k_vode, nr_vode, rho_src_vode, e_src_vode,&
                               rho_init_vode

    use eos_params_module
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_omb

    use cvode_interface
    use fnvector_serial
    use eos_module, only: vode_rtol, vode_atol_scaled
    use, intrinsic :: iso_c_binding

    implicit none

    real(rt), intent(in   ) :: dt
    real(rt), intent(in   ) :: rho_in, T_in, ne_in, e_in
    type(c_ptr), value :: cvmem
    type(c_ptr), value :: sunvec_y
    real(rt), intent(  out) ::  rho_out, T_out,ne_out,e_out
    real(rt),  intent(in   ) ::         rho_src, e_src

    real(rt) mean_rhob

    real(c_double), pointer, dimension(:) :: atol
    real(c_double) :: rtol
    real(c_double) :: time, tout
    integer(c_long), intent(in) :: neq
    real(c_double), pointer, intent(in) :: yvec(:)
    
    integer(c_int) :: ierr, maxord
    integer(c_long) :: mxsteps
    real(c_double) :: t_soln
    type(c_ptr) :: sunvec_atol
    logical, save :: firstCall = .true.

    allocate(atol(neq))

    t_vode   = t_in
    ne_vode  = ne_in
    rho_vode = rho_in
    rho_init_vode = rho_in
    nr_vode  = 0
    rho_src_vode = rho_src
    e_src_vode = e_src

    ! Initialize the integration time
    time = 0.d0
    
    ! We will integrate "e" and "rho" in time. 
    yvec(1) = e_in
    yvec(2) = rho_in

    ! Set the tolerances.  
    atol(1) = 1.d-4 * e_in
    atol(2) = 1.d-4 * rho_in
    rtol = 1.d-4

    maxord = 1
    mxsteps = 1

    sunvec_atol = n_vmake_serial(neq, atol)
    ierr = fcvodereinit(cvmem, time, sunvec_y)
    ierr = fcvodesvtolerances(cvmem, rtol, sunvec_atol)

    ierr = fcvodesetmaxord(cvmem, maxord) 
    ierr = fcvodesetmaxnumsteps(cvmem, mxsteps) 
    ierr = fcvodesetinitstep(cvmem, dt) 
    
    ierr = fcvode(cvmem, dt, sunvec_y, time, cv_normal)
    
    e_out  = yvec(1)
    rho_out = yvec(2)
!    rho_out = rho_init_vode + dt * rho_src_vode
    t_out  = t_vode
    ne_out = ne_vode

    if (ierr .ne. 0) then
       print *, 'istate = ', ierr, 'at (i,j,k) ',i_vode,j_vode,k_vode
       call amrex_error("ERROR in fcvode_wrapper_with_src: integration failed")
    endif

    call n_vdestroy_serial(sunvec_atol)

end subroutine fcvode_wrapper_with_source_single

    integer(c_int) function rhsfn_src(tn, sunvec_y, sunvec_f, user_data) &
           result(ierr) bind(c,name='RhsFn_src')

      use, intrinsic :: iso_c_binding
      use fnvector_serial
      use cvode_interface
      implicit none

      real(c_double), value :: tn
      type(c_ptr), value    :: sunvec_y, sunvec_f, user_data

      ! pointers to data in SUNDAILS vectors
      real(c_double), dimension(:), pointer :: yvec, fvec

      integer(c_long) :: neq
      real(c_double) :: energy(2)

      neq = int(2, c_long)

      ! get data arrays from SUNDIALS vectors
      call n_vgetdata_serial(sunvec_y, neq, yvec)
      call n_vgetdata_serial(sunvec_f, neq, fvec)

      call f_rhs_split(2, tn, yvec, energy, 0.0, 0)

      fvec = energy

      ierr = 0
    end function rhsfn_src
end module fcvode_extras_src

````

