
# File farkode\_extras.f90

[**File List**](files.md) **>** [**HeatCool**](dir_8c890215953ac09098af8cb94c8b9fc0.md) **>** [**farkode\_extras.f90**](farkode__extras_8f90.md)

[Go to the documentation of this file.](farkode__extras_8f90.md) 


````cpp
module farkode_extras

  implicit none

  contains

    subroutine farkode_wrapper(dt, rho_in, T_in, ne_in, e_in, neq, cvmem, &
                              sunvec_y, yvec, T_out, ne_out, e_out)

        use amrex_fort_module, only : rt => amrex_real
        use vode_aux_module, only: rho_vode, t_vode, ne_vode, z_vode
        use atomic_rates_module, only: this_z
        use arkode_interface
        use fnvector_serial
        use eos_module, only: vode_rtol, vode_atol_scaled
        use, intrinsic :: iso_c_binding

        implicit none

        real(rt), intent(in   ) :: dt
        real(rt), intent(in   ) :: rho_in, T_in, ne_in, e_in
        type(c_ptr), value :: cvmem
        type(c_ptr), value :: sunvec_y
        real(rt), intent(  out) ::         T_out,ne_out,e_out

        real(c_double) :: atol, rtol
        real(c_double) :: time, tout
        integer(c_long), intent(in) :: neq
        real(c_double), pointer, intent(in) :: yvec(:)

        integer(c_int) :: ierr

        real(c_double) :: t_soln

        t_vode   = t_in
        ne_vode  = ne_in
        rho_vode = rho_in

        ! Initialize the integration time
        time = 0.d0

        ! We will integrate "e" in time. 
        yvec(1) = e_in

        ! Set the tolerances.  
        atol = vode_atol_scaled * e_in
        rtol = vode_rtol

        ierr = farkodereinit(cvmem, c_null_ptr, c_funloc(rhsfnark), time, sunvec_y)
        ierr = farkodesstolerances(cvmem, rtol, atol)

        ierr = farkode(cvmem, dt, sunvec_y, time, ark_normal)

        e_out  = yvec(1)
        t_out  = t_vode
        ne_out = ne_vode

    end subroutine farkode_wrapper

    integer(c_int) function rhsfnark(tn, sunvec_y, sunvec_f, user_data) &
           result(ierr) bind(c,name='RhsFnArk')

      use, intrinsic :: iso_c_binding
      use fnvector_serial
      use cvode_interface
      implicit none

      real(c_double), value :: tn
      type(c_ptr), value    :: sunvec_y
      type(c_ptr), value    :: sunvec_f
      type(c_ptr), value    :: user_data

      ! pointers to data in SUNDAILS vectors
      real(c_double), pointer :: yvec(:)
      real(c_double), pointer :: fvec(:)

      real(c_double) :: energy

      integer(c_long), parameter :: neq = 1

      ! get data arrays from SUNDIALS vectors
      call n_vgetdata_serial(sunvec_y, neq, yvec)
      call n_vgetdata_serial(sunvec_f, neq, fvec)

      call f_rhs(1, tn, yvec(1), energy, 0.0, 0)

      fvec(1) = energy

      ierr = 0
    end function rhsfnark

end module farkode_extras
````

