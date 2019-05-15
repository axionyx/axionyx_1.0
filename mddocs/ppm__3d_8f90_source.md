
# File ppm\_3d.f90

[**File List**](files.md) **>** [**HydroFortran**](dir_1fab266cd447ad3f3624320661f845f1.md) **>** [**ppm\_3d.f90**](ppm__3d_8f90.md)

[Go to the documentation of this file.](ppm__3d_8f90.md) 


````cpp
module ppm_module

  implicit none

  private

  public ppm

contains
  !
  ! characteristics based on u
  !
  subroutine ppm(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                 u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                 flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                 Ip,Im, &
                 ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,a_old)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : ppm_type

    implicit none

    integer         , intent(in   ) ::   s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in   ) ::  qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer         , intent(in   ) ::   f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer         , intent(in   ) ::  ilo1,ilo2,ihi1,ihi2
 
    real(rt), intent(in   ) ::      s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    real(rt), intent(in   ) ::      u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    real(rt), intent(in   ) ::   cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt), intent(in   ) ::  flatn( f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    real(rt), intent(inout) :: Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    real(rt), intent(inout) :: Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    real(rt), intent(in   ) ::  dx,dy,dz,dt,a_old

    real(rt) :: dt_over_a
    integer          :: k3d,kc

    dt_over_a = dt / a_old
   
    if (ppm_type .eq. 1) then

        call ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       ip,im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

    else if (ppm_type .eq. 2) then

        call ppm_type2(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       ip,im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

    end if

  end subroutine ppm

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::
    
  subroutine ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

    use amrex_error_module
    use amrex_mempool_module, only: amrex_allocate, amrex_deallocate
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals
    use amrex_constants_module

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2

    real(rt), intent(in) :: s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    real(rt), intent(in) :: u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    real(rt), intent(in) :: cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt), intent(in) :: flatn(f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    real(rt), intent(out) :: Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    real(rt), intent(out) :: Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    ! Note that dt_over_a = dt / a_old
    real(rt), intent(in) :: dx,dy,dz,dt_over_a
    integer, intent(in)          :: k3d,kc

    ! local
    integer i,j,k

    real(rt) :: dxinv,dyinv,dzinv

    real(rt), pointer :: dsl(:), dsr(:), dsc(:)
    real(rt), pointer :: sigma(:), s6(:)

    ! s_{\ib,+}, s_{\ib,-}
    real(rt), pointer :: sp(:)
    real(rt), pointer :: sm(:)

    ! \delta s_{\ib}^{vL}
    real(rt), pointer :: dsvl(:,:)
    real(rt), pointer :: dsvlm(:,:)
    real(rt), pointer :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    real(rt), pointer :: sedge(:,:)
    real(rt), pointer :: sedgez(:,:,:)

    dxinv = 1.0d0/dx
    dyinv = 1.0d0/dy
    dzinv = 1.0d0/dz

    ! cell-centered indexing
    call amrex_allocate(sp,ilo1-1,ihi1+1)
    call amrex_allocate(sm,ilo1-1,ihi1+1)

    call amrex_allocate(sigma,ilo1-1,ihi1+1)
    call amrex_allocate(s6,ilo1-1,ihi1+1)

    if (ppm_type .ne. 1) &
         call amrex_error("Should have ppm_type = 1 in ppm_type1")

    if (s_l1 .gt. ilo1-3 .or. s_l2 .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_l1, s_l2
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if

    if (s_h1 .lt. ihi1+3 .or. s_h2 .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_h1, s_h2
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    call amrex_allocate(dsvl,ilo1-2,ihi1+2,ilo2-1,ihi2+1)

    ! edge-centered indexing for x-faces -- ppm_type = 1 only
    call amrex_allocate(sedge,ilo1-1,ihi1+2,ilo2-1,ihi2+1)

    ! cell-centered indexing
    call amrex_allocate(dsc,ilo1-2,ihi1+2)
    call amrex_allocate(dsl,ilo1-2,ihi1+2)
    call amrex_allocate(dsr,ilo1-2,ihi1+2)

    ! compute s at x-edges

    ! compute van Leer slopes in x-direction
    dsvl = zero
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+2
          dsc(i) = half * (s(i+1,j,k3d) - s(i-1,j,k3d))
          dsl(i) = two  * (s(i  ,j,k3d) - s(i-1,j,k3d))
          dsr(i) = two  * (s(i+1,j,k3d) - s(i  ,j,k3d))
          if (dsl(i)*dsr(i) .gt. zero) &
               dsvl(i,j) = sign(one,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! interpolate s to x-edges
       do i=ilo1-1,ihi1+2
          sedge(i,j) = half*(s(i,j,k3d)+s(i-1,j,k3d)) &
               - sixth*(dsvl(i,j)-dsvl(i-1,j))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i-1,j,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i-1,j,k3d)))
       end do

       ! copy sedge into sp and sm
       sp(ilo1-1:ihi1+1) = sedge(ilo1:ihi1+2,j)
       sm(ilo1-1:ihi1+1) = sedge(ilo1-1:ihi1+1,j)

       if (ppm_flatten_before_integrals == 1) then
          ! Flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! Modify using quadratic limiters -- note this version of the limiting comes
       ! from Colella and Sekora (2008), not the original PPM paper.
       do i = ilo1-1, ihi1+1
          if ((sp(i)-s(i,j,k3d))*(s(i,j,k3d)-sm(i)) .le. zero) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(sp(i)-s(i,j,k3d)) .ge. two*abs(sm(i)-s(i,j,k3d))) then
             sp(i) = three*s(i,j,k3d) - two*sm(i)
          else if (abs(sm(i)-s(i,j,k3d)) .ge. two*abs(sp(i)-s(i,j,k3d))) then
             sm(i) = three*s(i,j,k3d) - two*sp(i)
          end if
       end do

       if (ppm_flatten_before_integrals == 2) then
          ! Flatten the parabola AFTER doing the monotonization --
          ! this is the method that Miller & Colella do
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! compute x-component of Ip and Im
       s6 = six*s(ilo1-1:ihi1+1,j,k3d) - three*(sm+sp)

       ! Ip/m is the integral under the parabola for the extent
       ! that a wave can travel over a timestep
       !
       ! Ip integrates to the right edge of a cell
       ! Im integrates to the left edge of a cell

       ! u-c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,1)-cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dxinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)-cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,1,1) = sp(i)
          else
             ip(i,j,kc,1,1) = sp(i) - &
               half*sigma(i)*(sp(i)-sm(i)-(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)-cspd(i,j,k3d) >= zero) then
             im(i,j,kc,1,1) = sm(i)
          else
             im(i,j,kc,1,1) = sm(i) + &
               half*sigma(i)*(sp(i)-sm(i)+(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       ! u wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,1))*dt_over_a*dxinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1) <= zero) then
             ip(i,j,kc,1,2) = sp(i)
          else
             ip(i,j,kc,1,2) = sp(i) - &
               half*sigma(i)*(sp(i)-sm(i)-(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1) >= zero) then
             im(i,j,kc,1,2) = sm(i)
          else
             im(i,j,kc,1,2) = sm(i) + &
               half*sigma(i)*(sp(i)-sm(i)+(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       ! u+c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,1)+cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dxinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)+cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,1,3) = sp(i)
          else
             ip(i,j,kc,1,3) = sp(i) - &
               half*sigma(i)*(sp(i)-sm(i)-(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)+cspd(i,j,k3d) >= zero) then
             im(i,j,kc,1,3) = sm(i)
          else
             im(i,j,kc,1,3) = sm(i) + &
               half*sigma(i)*(sp(i)-sm(i)+(one-two3rd*sigma(i))*s6(i))
          endif
       end do

    end do

    call amrex_deallocate(dsc)
    call amrex_deallocate(dsl)
    call amrex_deallocate(dsr)
    call amrex_deallocate(sedge)
    call amrex_deallocate(dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    call amrex_allocate( dsvl,ilo1-1,ihi1+1,ilo2-2,ihi2+2)

    ! edge-centered indexing for y-faces
    call amrex_allocate(sedge,ilo1-1,ihi1+1,ilo2-1,ihi2+2)

    ! cell-centered indexing
    call amrex_allocate(dsc,ilo1-1,ihi1+1)
    call amrex_allocate(dsl,ilo1-1,ihi1+1)
    call amrex_allocate(dsr,ilo1-1,ihi1+1)

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    dsvl = zero
    do j=ilo2-2,ihi2+2
       do i=ilo1-1,ihi1+1
          dsc(i) = half * (s(i,j+1,k3d) - s(i,j-1,k3d))
          dsl(i) = two  * (s(i,j  ,k3d) - s(i,j-1,k3d))
          dsr(i) = two  * (s(i,j+1,k3d) - s(i,j  ,k3d))
          if (dsl(i)*dsr(i) .gt. zero) &
               dsvl(i,j) = sign(one,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do
    end do

    ! interpolate s to y-edges
    do j=ilo2-1,ihi2+2
       do i=ilo1-1,ihi1+1
          sedge(i,j) = half*(s(i,j,k3d)+s(i,j-1,k3d)) &
               - sixth*(dsvl(i,j)-dsvl(i,j-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i,j-1,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i,j-1,k3d)))
       end do
    end do

    do j=ilo2-1,ihi2+1

       ! copy sedge into sp and sm
       sp = sedge(ilo1-1:ihi1+1,j+1)
       sm = sedge(ilo1-1:ihi1+1,j  )

       if (ppm_flatten_before_integrals == 1) then
          ! Flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! Modify using quadratic limiters
       do i = ilo1-1, ihi1+1
          if ((sp(i)-s(i,j,k3d))*(s(i,j,k3d)-sm(i)) .le. zero) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(sp(i)-s(i,j,k3d)) .ge. two*abs(sm(i)-s(i,j,k3d))) then
             sp(i) = three*s(i,j,k3d) - two*sm(i)
          else if (abs(sm(i)-s(i,j,k3d)) .ge. two*abs(sp(i)-s(i,j,k3d))) then
             sm(i) = three*s(i,j,k3d) - two*sp(i)
          end if
       end do

       if (ppm_flatten_before_integrals == 2) then
          ! Flatten the parabola AFTER doing the monotonization --
          ! this is the method that Miller & Colella do
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! compute y-component of Ip and Im
       s6 = six*s(ilo1-1:ihi1+1,j,k3d) - three*(sm+sp)

       ! v-c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,2)-cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dyinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)-cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,2,1) = sp(i)
          else
             ip(i,j,kc,2,1) = sp(i) - &
               half*sigma(i)*(sp(i)-sm(i)-(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)-cspd(i,j,k3d) >= zero) then
             im(i,j,kc,2,1) = sm(i)
          else
             im(i,j,kc,2,1) = sm(i) + &
               half*sigma(i)*(sp(i)-sm(i)+(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       ! v wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,2))*dt_over_a*dyinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2) <= zero) then
             ip(i,j,kc,2,2) = sp(i)
          else
             ip(i,j,kc,2,2) = sp(i) - &
               half*sigma(i)*(sp(i)-sm(i)-(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2) >= zero) then
             im(i,j,kc,2,2) = sm(i)
          else
             im(i,j,kc,2,2) = sm(i) + &
               half*sigma(i)*(sp(i)-sm(i)+(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       ! v+c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,2)+cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dyinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)+cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,2,3) = sp(i)
          else
             ip(i,j,kc,2,3) = sp(i) - &
               half*sigma(i)*(sp(i)-sm(i)-(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)+cspd(i,j,k3d) >= zero) then
             im(i,j,kc,2,3) = sm(i)
          else
             im(i,j,kc,2,3) = sm(i) + &
               half*sigma(i)*(sp(i)-sm(i)+(one-two3rd*sigma(i))*s6(i))
          endif
       end do

    end do

    call amrex_deallocate(dsc)
    call amrex_deallocate(dsl)
    call amrex_deallocate(dsr)
    call amrex_deallocate(dsvl)
    call amrex_deallocate(sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    call amrex_allocate( dsvl,ilo1-1,ihi1+1,ilo2-1,ihi2+1)
    call amrex_allocate(dsvlm,ilo1-1,ihi1+1,ilo2-1,ihi2+1)
    call amrex_allocate(dsvlp,ilo1-1,ihi1+1,ilo2-1,ihi2+1)

    ! cell-centered indexing
    call amrex_allocate(dsc,ilo1-1,ihi1+1)
    call amrex_allocate(dsl,ilo1-1,ihi1+1)
    call amrex_allocate(dsr,ilo1-1,ihi1+1)

    call amrex_allocate(sedgez,ilo1-1,ihi1+1,ilo2-2,ihi2+3,k3d-1,k3d+2)

    ! compute s at z-edges

    ! compute van Leer slopes in z-direction
    dsvl  = zero
    dsvlm = zero
    dsvlp = zero

    do j=ilo2-1,ihi2+1

       ! compute on slab below
       k = k3d-1
       dsc = half * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k-1))
       dsl = two  * (s(ilo1-1:ihi1+1,j,k  ) - s(ilo1-1:ihi1+1,j,k-1))
       dsr = two  * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k  ))
       do i = ilo1-1, ihi1+1
          if (dsl(i)*dsr(i) .gt. zero) &
               dsvlm(i,j) = sign(one,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! compute on slab above
       k = k3d+1
       dsc = half * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k-1))
       dsl = two  * (s(ilo1-1:ihi1+1,j,k  ) - s(ilo1-1:ihi1+1,j,k-1))
       dsr = two  * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k  ))
       do i = ilo1-1, ihi1+1
          if (dsl(i)*dsr(i) .gt. zero) &
               dsvlp(i,j) = sign(one,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! compute on current slab
       k = k3d
       dsc = half * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k-1))
       dsl = two  * (s(ilo1-1:ihi1+1,j,k  ) - s(ilo1-1:ihi1+1,j,k-1))
       dsr = two  * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k  ))
       do i = ilo1-1, ihi1+1
          if (dsl(i)*dsr(i) .gt. zero) &
               dsvl(i,j) = sign(one,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! interpolate to lo face
       k = k3d
       sm = half*(s(ilo1-1:ihi1+1,j,k)+s(ilo1-1:ihi1+1,j,k-1)) - sixth*(dsvl(ilo1-1:ihi1+1,j)-dsvlm(ilo1-1:ihi1+1,j))
       ! make sure sedge lies in between adjacent cell-centered values
       sm = max(sm,min(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))
       sm = min(sm,max(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))

       ! interpolate to hi face
       k = k3d+1
       sp = half*(s(ilo1-1:ihi1+1,j,k)+s(ilo1-1:ihi1+1,j,k-1)) - sixth*(dsvlp(ilo1-1:ihi1+1,j)-dsvl(ilo1-1:ihi1+1,j))
       ! make sure sedge lies in between adjacent cell-centered values
       sp = max(sp,min(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))
       sp = min(sp,max(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))

       if (ppm_flatten_before_integrals == 1) then
          ! flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! modify using quadratic limiters
       do i = ilo1-1, ihi1+1
          if ((sp(i)-s(i,j,k3d))*(s(i,j,k3d)-sm(i)) .le. zero) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(sp(i)-s(i,j,k3d)) .ge. two*abs(sm(i)-s(i,j,k3d))) then
             sp(i) = three*s(i,j,k3d) - two*sm(i)
          else if (abs(sm(i)-s(i,j,k3d)) .ge. two*abs(sp(i)-s(i,j,k3d))) then
             sm(i) = three*s(i,j,k3d) - two*sp(i)
          end if
       end do

       if (ppm_flatten_before_integrals == 2) then
          ! flatten the parabola AFTER doing the monotonization --
          ! this is the method that Miller & Colella do
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (one-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! compute z-component of Ip and Im
       s6 = six*s(ilo1-1:ihi1+1,j,k3d) - three*(sm+sp)

       ! w-c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,3)-cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dzinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)-cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,3,1) = sp(i)
          else
             ip(i,j,kc,3,1) = sp(i) - &
               half*sigma(i)*(sp(i)-sm(i)-(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)-cspd(i,j,k3d) >= zero) then
             im(i,j,kc,3,1) = sm(i)
          else
             im(i,j,kc,3,1) = sm(i) + &
               half*sigma(i)*(sp(i)-sm(i)+(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       ! w wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,3))*dt_over_a*dzinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3) <= zero) then
             ip(i,j,kc,3,2) = sp(i)
          else
             ip(i,j,kc,3,2) = sp(i) - &
               half*sigma(i)*(sp(i)-sm(i)-(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3) >= zero) then
             im(i,j,kc,3,2) = sm(i)
          else
             im(i,j,kc,3,2) = sm(i) + &
               half*sigma(i)*(sp(i)-sm(i)+(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       ! w+c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,3)+cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dzinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)+cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,3,3) = sp(i)
          else
             ip(i,j,kc,3,3) = sp(i) - &
               half*sigma(i)*(sp(i)-sm(i)-(one-two3rd*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)+cspd(i,j,k3d) >= zero) then
             im(i,j,kc,3,3) = sm(i)
          else
             im(i,j,kc,3,3) = sm(i) + &
               half*sigma(i)*(sp(i)-sm(i)+(one-two3rd*sigma(i))*s6(i))
          endif
       end do

    end do

    call amrex_deallocate(dsc)
    call amrex_deallocate(dsl)
    call amrex_deallocate(dsr)
    call amrex_deallocate(dsvl)
    call amrex_deallocate(dsvlm)
    call amrex_deallocate(dsvlp)
    call amrex_deallocate(sp)
    call amrex_deallocate(sm)
    call amrex_deallocate(sedgez)
    call amrex_deallocate(sigma)
    call amrex_deallocate(s6)

  end subroutine ppm_type1

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine ppm_type2(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals
    use amrex_constants_module

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2

    real(rt)    s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    real(rt)    u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    real(rt) cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt) flatn( f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    real(rt) Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    real(rt) Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    ! Note that dt_over_a = dt / a_old
    real(rt) dx,dy,dz,dt_over_a
    real(rt) dxinv,dyinv,dzinv
    integer          k3d,kc

    ! local
    integer i,j,k
    logical extremum, bigp, bigm

    real(rt) D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    real(rt) sgn, sigma, s6
    real(rt) dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
    real(rt) dachkm, dachkp
    real(rt) amax, delam, delap

    ! s_{\ib,+}, s_{\ib,-}
    real(rt), allocatable :: sp(:,:)
    real(rt), allocatable :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    real(rt), allocatable :: dsvl(:,:)
    real(rt), allocatable :: dsvlm(:,:)
    real(rt), allocatable :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    real(rt), allocatable :: sedge(:,:)
    real(rt), allocatable :: sedgez(:,:,:)

    ! constant used in Colella 2008
    real(rt), parameter :: C = 1.25d0

    dxinv = 1.0d0/dx
    dyinv = 1.0d0/dy
    dzinv = 1.0d0/dz

    ! cell-centered indexing
    allocate(sp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(sm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    if (ppm_type .ne. 2) &
         call amrex_error("Should have ppm_type = 2 in ppm_type2")

    if (s_l1 .gt. ilo1-3 .or. s_l2 .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_l1, s_l2
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call amrex_error("Need more ghost cells on array in ppm_type2")
    end if

    if (s_h1 .lt. ihi1+3 .or. s_h2 .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_h1, s_h2
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call amrex_error("Need more ghost cells on array in ppm_type2")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(ilo1-2:ihi1+2,ilo2-1:ihi2+1))

    ! edge-centered indexing for x-faces
    allocate(sedge(ilo1-2:ihi1+3,ilo2-1:ihi2+1))

    ! compute s at x-edges

    ! interpolate s to x-edges
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+3
          sedge(i,j) = seven12th*(s(i-1,j,k3d)+s(i  ,j,k3d)) &
               - twelfth*(s(i-2,j,k3d)+s(i+1,j,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i-1,j,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. zero) then
             d2  = three*(s(i-1,j,k3d)-two*sedge(i,j)+s(i,j,k3d))
             d2l = s(i-2,j,k3d)-two*s(i-1,j,k3d)+s(i,j,k3d)
             d2r = s(i-1,j,k3d)-two*s(i,j,k3d)+s(i+1,j,k3d)
             sgn = sign(one,d2)
             d2lim = sgn*max(min(c*sgn*d2l,c*sgn*d2r,sgn*d2),zero)
             sedge(i,j) = half*(s(i-1,j,k3d)+s(i,j,k3d)) - sixth*d2lim
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i+1,j)-s(i,j,k3d)
          alpham   = sedge(i  ,j)-s(i,j,k3d)
          bigp     = abs(alphap).gt.two*abs(alpham)
          bigm     = abs(alpham).gt.two*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. zero) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i-1,j)
             dafacep   = sedge(i+2,j) - sedge(i+1,j)
             dabarm    = s(i,j,k3d) - s(i-1,j,k3d)
             dabarp    = s(i+1,j,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. zero)
          end if

          if (extremum) then
             d2     = six*(alpham + alphap)
             d2l    = s(i-2,j,k3d)-two*s(i-1,j,k3d)+s(i,j,k3d)
             d2r    = s(i,j,k3d)-two*s(i+1,j,k3d)+s(i+2,j,k3d)
             d2c    = s(i-1,j,k3d)-two*s(i,j,k3d)+s(i+1,j,k3d)
             sgn    = sign(one,d2)
             d2lim  = max(min(sgn*d2,c*sgn*d2l,c*sgn*d2r,c*sgn*d2c),zero)
             alpham = alpham*d2lim/max(abs(d2),1.d-10)
             alphap = alphap*d2lim/max(abs(d2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(one,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-two*delam - two*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -two*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(one,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-two*delap - two*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -two*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap

          if (ppm_flatten_before_integrals > 0) then
             ! flatten the parabola AFTER doing the monotonization
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (one-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (one-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          !
          ! Compute x-component of Ip and Im.
          !
          s6    = six*s(i,j,k3d) - three*(sm(i,j)+sp(i,j))

          ! u-c wave
          sigma = abs(u(i,j,k3d,1)-cspd(i,j,k3d))*dt_over_a*dxinv

          if (u(i,j,k3d,1)-cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,1,1) = sp(i,j)
          else
             ip(i,j,kc,1,1) = sp(i,j) - &
                  half*sigma*(sp(i,j)-sm(i,j)-(one-two3rd*sigma)*s6)
          endif

          if (u(i,j,k3d,1)-cspd(i,j,k3d) >= zero) then
             im(i,j,kc,1,1) = sm(i,j)
          else
             im(i,j,kc,1,1) = sm(i,j) + &
                  half*sigma*(sp(i,j)-sm(i,j)+(one-two3rd*sigma)*s6)
          endif

          ! u wave
          sigma = abs(u(i,j,k3d,1))*dt_over_a*dxinv

          if (u(i,j,k3d,1) <= zero) then
             ip(i,j,kc,1,2) = sp(i,j)
          else
             ip(i,j,kc,1,2) = sp(i,j) - &
                  half*sigma*(sp(i,j)-sm(i,j)-(one-two3rd*sigma)*s6)
          endif

          if (u(i,j,k3d,1) >= zero) then
             im(i,j,kc,1,2) = sm(i,j)
          else
             im(i,j,kc,1,2) = sm(i,j) + &
                  half*sigma*(sp(i,j)-sm(i,j)+(one-two3rd*sigma)*s6)
          endif

          ! u+c wave
          sigma = abs(u(i,j,k3d,1)+cspd(i,j,k3d))*dt_over_a*dxinv

          if (u(i,j,k3d,1)+cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,1,3) = sp(i,j) 
          else
             ip(i,j,kc,1,3) = sp(i,j) - &
                  half*sigma*(sp(i,j)-sm(i,j)-(one-two3rd*sigma)*s6)
          endif

          if (u(i,j,k3d,1)+cspd(i,j,k3d) >= zero) then
             im(i,j,kc,1,3) = sm(i,j) 
          else
             im(i,j,kc,1,3) = sm(i,j) + &
                  half*sigma*(sp(i,j)-sm(i,j)+(one-two3rd*sigma)*s6)
          endif

       end do
    end do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-2:ihi2+2))

    ! edge-centered indexing for y-faces
    allocate(sedge(ilo1-1:ihi1+1,ilo2-2:ihi2+3))

    ! compute s at y-edges

    ! interpolate s to y-edges
    do j=ilo2-2,ihi2+3
       do i=ilo1-1,ihi1+1
          sedge(i,j) = seven12th*(s(i,j-1,k3d)+s(i,j,k3d)) &
               - twelfth*(s(i,j-2,k3d)+s(i,j+1,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i,j-1,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. zero) then
             d2  = three*(s(i,j-1,k3d)-two*sedge(i,j)+s(i,j,k3d))
             d2l = s(i,j-2,k3d)-two*s(i,j-1,k3d)+s(i,j,k3d)
             d2r = s(i,j-1,k3d)-two*s(i,j,k3d)+s(i,j+1,k3d)
             sgn = sign(one,d2)
             d2lim = sgn*max(min(c*sgn*d2l,c*sgn*d2r,sgn*d2),zero)
             sedge(i,j) = half*(s(i,j-1,k3d)+s(i,j,k3d)) - sixth*d2lim
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i,j+1)-s(i,j,k3d)
          alpham   = sedge(i,j  )-s(i,j,k3d)
          bigp     = abs(alphap).gt.two*abs(alpham)
          bigm     = abs(alpham).gt.two*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. zero) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i,j-1)
             dafacep   = sedge(i,j+2) - sedge(i,j+1)
             dabarm    = s(i,j,k3d) - s(i,j-1,k3d)
             dabarp    = s(i,j+1,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. zero)
          end if

          if (extremum) then
             d2     = six*(alpham + alphap)
             d2l    = s(i,j-2,k3d)-two*s(i,j-1,k3d)+s(i,j,k3d)
             d2r    = s(i,j,k3d)-two*s(i,j+1,k3d)+s(i,j+2,k3d)
             d2c    = s(i,j-1,k3d)-two*s(i,j,k3d)+s(i,j+1,k3d)
             sgn    = sign(one,d2)
             d2lim  = max(min(sgn*d2,c*sgn*d2l,c*sgn*d2r,c*sgn*d2c),zero)
             alpham = alpham*d2lim/max(abs(d2),1.d-10)
             alphap = alphap*d2lim/max(abs(d2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(one,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j-1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-two*delam - two*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -two*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(one,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j+1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-two*delap - two*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -two*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap

          if (ppm_flatten_before_integrals > 0) then
             ! flatten the parabola AFTER doing the monotonization
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (one-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (one-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          !
          ! Compute y-component of Ip and Im.
          !
          s6    = six*s(i,j,k3d) - three*(sm(i,j)+sp(i,j))

          ! v-c wave
          sigma = abs(u(i,j,k3d,2)-cspd(i,j,k3d))*dt_over_a*dyinv

          if (u(i,j,k3d,2)-cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,2,1) = sp(i,j) 
          else
             ip(i,j,kc,2,1) = sp(i,j) - &
                  half*sigma*(sp(i,j)-sm(i,j)-(one-two3rd*sigma)*s6)
          endif

          if (u(i,j,k3d,2)-cspd(i,j,k3d) >= zero) then
             im(i,j,kc,2,1) = sm(i,j) 
          else
             im(i,j,kc,2,1) = sm(i,j) + &
                  half*sigma*(sp(i,j)-sm(i,j)+(one-two3rd*sigma)*s6)
          endif

          ! v wave
          sigma = abs(u(i,j,k3d,2))*dt_over_a*dyinv

          if (u(i,j,k3d,2) <= zero) then
             ip(i,j,kc,2,2) = sp(i,j) 
          else
             ip(i,j,kc,2,2) = sp(i,j) - &
                  half*sigma*(sp(i,j)-sm(i,j)-(one-two3rd*sigma)*s6)
          endif

          if (u(i,j,k3d,2) >= zero) then
             im(i,j,kc,2,2) = sm(i,j) 
          else
             im(i,j,kc,2,2) = sm(i,j) + &
                  half*sigma*(sp(i,j)-sm(i,j)+(one-two3rd*sigma)*s6)
          endif

          ! v+c wave
          sigma = abs(u(i,j,k3d,2)+cspd(i,j,k3d))*dt_over_a*dyinv

          if (u(i,j,k3d,2)+cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,2,3) = sp(i,j) 
          else
             ip(i,j,kc,2,3) = sp(i,j) - &
                  half*sigma*(sp(i,j)-sm(i,j)-(one-two3rd*sigma)*s6)
          endif

          if (u(i,j,k3d,2)+cspd(i,j,k3d) >= zero) then
             im(i,j,kc,2,3) = sm(i,j) 
          else
             im(i,j,kc,2,3) = sm(i,j) + &
                  half*sigma*(sp(i,j)-sm(i,j)+(one-two3rd*sigma)*s6)
          endif
       end do
    end do

    deallocate(dsvl,sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    allocate(sedgez(ilo1-1:ihi1+1,ilo2-2:ihi2+3,k3d-1:k3d+2))

    ! compute s at z-edges

    ! interpolate s to z-edges
    do k=k3d-1,k3d+2
       do j=ilo2-1,ihi2+1
          do i=ilo1-1,ihi1+1
             sedgez(i,j,k) = seven12th*(s(i,j,k-1)+s(i,j,k)) &
                  - twelfth*(s(i,j,k-2)+s(i,j,k+1))
             !
             ! limit sedgez
             !
             if ((sedgez(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedgez(i,j,k)) .lt. zero) then
                d2  = three*(s(i,j,k-1)-two*sedgez(i,j,k)+s(i,j,k))
                d2l = s(i,j,k-2)-two*s(i,j,k-1)+s(i,j,k)
                d2r = s(i,j,k-1)-two*s(i,j,k)+s(i,j,k+1)
                sgn = sign(one,d2)
                d2lim = sgn*max(min(c*sgn*d2l,c*sgn*d2r,sgn*d2),zero)
                sedgez(i,j,k) = half*(s(i,j,k-1)+s(i,j,k)) - sixth*d2lim
             end if
          end do
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    k = k3d
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedgez(i,j,k+1)-s(i,j,k)
          alpham   = sedgez(i,j,k  )-s(i,j,k)
          bigp     = abs(alphap).gt.two*abs(alpham)
          bigm     = abs(alpham).gt.two*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. zero) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedgez(i,j,k) - sedgez(i,j,k-1)
             dafacep   = sedgez(i,j,k+2) - sedgez(i,j,k+1)
             dabarm    = s(i,j,k) - s(i,j,k-1)
             dabarp    = s(i,j,k+1) - s(i,j,k)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. zero)
          end if

          if (extremum) then
             d2     = six*(alpham + alphap)
             d2l    = s(i,j,k-2)-two*s(i,j,k-1)+s(i,j,k)
             d2r    = s(i,j,k)-two*s(i,j,k+1)+s(i,j,k+2)
             d2c    = s(i,j,k-1)-two*s(i,j,k)+s(i,j,k+1)
             sgn    = sign(one,d2)
             d2lim  = max(min(sgn*d2,c*sgn*d2l,c*sgn*d2r,c*sgn*d2c),zero)
             alpham = alpham*d2lim/max(abs(d2),1.d-10)
             alphap = alphap*d2lim/max(abs(d2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(one,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j,k-1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-two*delam - two*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -two*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(one,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j,k+1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-two*delap - two*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -two*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k) + alpham
          sp(i,j) = s(i,j,k) + alphap

          if (ppm_flatten_before_integrals > 0) then
             ! flatten the parabola AFTER doing the monotonization (note k = k3d here)
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (one-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (one-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          !
          ! Compute z-component of Ip and Im.
          !
          s6    = six*s(i,j,k3d) - three*(sm(i,j)+sp(i,j))
          
          ! w-c wave
          sigma = abs(u(i,j,k3d,3)-cspd(i,j,k3d))*dt_over_a*dzinv
          
          if (u(i,j,k3d,3)-cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,3,1) = sp(i,j) 
          else
             ip(i,j,kc,3,1) = sp(i,j) - &
                  half*sigma*(sp(i,j)-sm(i,j)-(one-two3rd*sigma)*s6)
          endif

          if (u(i,j,k3d,3)-cspd(i,j,k3d) >= zero) then
             im(i,j,kc,3,1) = sm(i,j) 
          else
             im(i,j,kc,3,1) = sm(i,j) + &
                  half*sigma*(sp(i,j)-sm(i,j)+(one-two3rd*sigma)*s6)
          endif

          ! w wave
          sigma = abs(u(i,j,k3d,3))*dt_over_a*dzinv

          if (u(i,j,k3d,3) <= zero) then
             ip(i,j,kc,3,2) = sp(i,j)
          else
             ip(i,j,kc,3,2) = sp(i,j) - &
                  half*sigma*(sp(i,j)-sm(i,j)-(one-two3rd*sigma)*s6)
          endif

          if (u(i,j,k3d,3) >= zero) then
             im(i,j,kc,3,2) = sm(i,j)
          else
             im(i,j,kc,3,2) = sm(i,j) + &
                  half*sigma*(sp(i,j)-sm(i,j)+(one-two3rd*sigma)*s6)
          endif

          ! w+c wave
          sigma = abs(u(i,j,k3d,3)+cspd(i,j,k3d))*dt_over_a*dzinv

          if (u(i,j,k3d,3)+cspd(i,j,k3d) <= zero) then
             ip(i,j,kc,3,3) = sp(i,j) 
          else
             ip(i,j,kc,3,3) = sp(i,j) - &
                  half*sigma*(sp(i,j)-sm(i,j)-(one-two3rd*sigma)*s6)
          endif

          if (u(i,j,k3d,3)+cspd(i,j,k3d) >= zero) then
             im(i,j,kc,3,3) = sm(i,j) 
          else
             im(i,j,kc,3,3) = sm(i,j) + &
                  half*sigma*(sp(i,j)-sm(i,j)+(one-two3rd*sigma)*s6)
          endif

       end do
    end do

    deallocate(dsvl,dsvlm,dsvlp,sp,sm,sedgez)

  end subroutine ppm_type2

end module ppm_module

````

