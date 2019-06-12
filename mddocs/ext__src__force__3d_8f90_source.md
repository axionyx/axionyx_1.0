
# File ext\_src\_force\_3d.f90

[**File List**](files.md) **>** [**Forcing**](dir_45682215f16eaf57f766b3c547de68bc.md) **>** [**ext\_src\_force\_3d.f90**](ext__src__force__3d_8f90.md)

[Go to the documentation of this file.](ext__src__force__3d_8f90.md) 


````cpp
subroutine ext_src_force(lo, hi, old_state, os_l1, os_l2, os_l3, os_h1, os_h2, os_h3, &
                                 new_state, ns_l1, ns_l2, ns_l3, ns_h1, ns_h2, ns_h3, &
                                 old_diag , od_l1, od_l2, od_l3, od_h1, od_h2, od_h3, &
                                 new_diag , nd_l1, nd_l2, nd_l3, nd_h1, nd_h2, nd_h3, &
                         src, src_l1, &
                         src_l2, src_l3, src_h1, src_h2, src_h3, problo, dx, time, z, dt)
!
!   Calculates the sources to be added later on.
!
!   Parameters
!   ----------
!   lo : double array (3)
!       The low corner of the current box.
!   hi : double array (3)
!       The high corner of the current box.
!   old_state_* : double arrays
!       The state vars
!   new_state_* : double arrays
!       The state vars
!   src_* : doubles arrays
!       The source terms to be added to state (iterative approx.)
!   problo : double array (3)
!       The low corner of the entire domain
!   dx : double array (3)
!       The cell size of this level.
!   time : double
!       The current time, in Mpc km^-1 s ~ 10^12 yr.
!      z : double
!       The current z = 1 / a - 1
!   dt : double
!       time step size, in Mpc km^-1 s ~ 10^12 yr.
!
!   Returns
!   -------
!   src : double array (dims) @todo
!       @todo
!
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : nvar, umx, umy, umz, ueden, ueint
    use fundamental_constants_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: os_l1, os_l2, os_l3
    integer, intent(in) :: os_h1, os_h2, os_h3
    integer, intent(in) :: ns_l1, ns_l2, ns_l3
    integer, intent(in) :: ns_h1, ns_h2, ns_h3
    integer, intent(in) :: od_l1, od_l2, od_l3
    integer, intent(in) :: od_h1, od_h2, od_h3
    integer, intent(in) :: nd_l1, nd_l2, nd_l3
    integer, intent(in) :: nd_h1, nd_h2, nd_h3
    integer, intent(in) :: src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
    real(rt), intent(in) :: old_state(os_l1:os_h1, &
                                              os_l2:os_h2, &
                                              os_l3:os_h3, NVAR)
    real(rt), intent(in) :: new_state(ns_l1:ns_h1, &
                                              ns_l2:ns_h2, &
                                              ns_l3:ns_h3, NVAR)
    real(rt), intent(in) :: old_diag (od_l1:od_h1, &
                                              od_l2:od_h2, &
                                              od_l3:od_h3, 2)
    real(rt), intent(in) :: new_diag (nd_l1:nd_h1, &
                                              nd_l2:nd_h2, &
                                              nd_l3:nd_h3, 2)
    real(rt), intent(in) :: problo(3), dx(3), z, dt, time

    real(rt), intent(out) :: src(src_l1:src_h1, src_l2:src_h2, &
                                         src_l3:src_h3, NVAR)

    real(rt), allocatable :: tmp_state(:,:,:,:)

    integer  :: i, j, k
    real(rt) :: a, half_dt

    ! Make a copy of the state so we can evolve it then throw it away
    allocate(tmp_state(ns_l1:ns_h1,ns_l2:ns_h2,ns_l3:ns_h3,nvar))
    tmp_state(:,:,:,:) = new_state(:,:,:,:)

    a = 1.d0 / (1.d0+z)

    ! Note that when we call this routine to compute the "old" source,
    !      both "old_state" and "new_state" are acutally the "old" state.
    ! When we call this routine to compute the "new" source,
    !      both "old_state" is in fact the "old" state and
    !           "new_state" is in fact the "new" state

    half_dt = 0.5d0 * dt

    call integrate_state_force(lo,hi,tmp_state,ns_l1,ns_l2,ns_l3,ns_h1,ns_h2,ns_h3, &
                                     new_diag ,nd_l1,nd_l2,nd_l3,nd_h1,nd_h2,nd_h3, &
                               dx,time,a,half_dt)

    ! Recall that this routine is called from a tiled MFIter 
    !   For old source: lo(:), hi(:) are the bounds of the growntilebox(src.nGrow9))
    !   For new source: lo(:), hi(:) are the bounds of the      tilebox, e.g. valid region only
    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               src(i,j,k,umx)   = (tmp_state(i,j,k,umx)   - new_state(i,j,k,umx)) * a / half_dt
               src(i,j,k,umy)   = (tmp_state(i,j,k,umy)   - new_state(i,j,k,umy)) * a / half_dt
               src(i,j,k,umz)   = (tmp_state(i,j,k,umz)   - new_state(i,j,k,umz)) * a / half_dt
               src(i,j,k,ueint) = (tmp_state(i,j,k,ueint) - new_state(i,j,k,ueint)) * a / half_dt
               src(i,j,k,ueden) = (tmp_state(i,j,k,ueden) - new_state(i,j,k,ueden)) * a / half_dt
           end do
        end do
    end do

end subroutine ext_src_force
````

