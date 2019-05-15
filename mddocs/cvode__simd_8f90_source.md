
# File cvode\_simd.f90

[**File List**](files.md) **>** [**Initialization**](dir_71a4420ed1f8982e7234eb6a0b7e6d5d.md) **>** [**cvode\_simd.f90**](cvode__simd_8f90.md)

[Go to the documentation of this file.](cvode__simd_8f90.md) 


````cpp
subroutine set_simd (simd_width_in) bind(C, name='set_simd')

   use misc_params, only: simd_width
   implicit none

   integer, intent(in) :: simd_width_in

   simd_width = simd_width_in

end subroutine set_simd

subroutine fort_alloc_simd_vec() bind(C, name='fort_alloc_simd_vec')
  use misc_params, only: simd_width
  use vode_aux_module, only: t_vode_vec, ne_vode_vec, rho_vode_vec
  use amrex_error_module, only: amrex_abort
  implicit none

  !$omp parallel
  if (allocated(t_vode_vec) .or. allocated(ne_vode_vec) .or. allocated(rho_vode_vec)) then
    !$omp single
    call amrex_abort("Why are VODE SIMD vectors already allocated??")
    !$omp end single
  end if

  allocate(t_vode_vec(simd_width), ne_vode_vec(simd_width), rho_vode_vec(simd_width))
  !$omp end parallel
end subroutine fort_alloc_simd_vec


subroutine fort_dealloc_simd_vec() bind(C, name='fort_dealloc_simd_vec')
  use vode_aux_module, only: t_vode_vec, ne_vode_vec, rho_vode_vec
  use amrex_error_module, only: amrex_abort
  implicit none

  !$omp parallel
  if (.not. (allocated(t_vode_vec) .and. allocated(ne_vode_vec) .and. allocated(rho_vode_vec))) then
    !$omp single
    call amrex_abort("Why are VODE SIMD vectors already deallocated??")
    !$omp end single
  end if

  deallocate(t_vode_vec, ne_vode_vec, rho_vode_vec)
  !$omp end parallel
end subroutine fort_dealloc_simd_vec
````

