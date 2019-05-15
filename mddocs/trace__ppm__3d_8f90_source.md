
# File trace\_ppm\_3d.f90

[**File List**](files.md) **>** [**HydroFortran**](dir_1fab266cd447ad3f3624320661f845f1.md) **>** [**trace\_ppm\_3d.f90**](trace__ppm__3d_8f90.md)

[Go to the documentation of this file.](trace__ppm__3d_8f90.md) 


````cpp
module trace_ppm_module

  implicit none

  private

  public tracexy_ppm, tracez_ppm

contains

    subroutine tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           Ip,Im,Ip_g,Im_g, &
                           qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                           srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                           ilo1,ilo2,ihi1,ihi2,dt,a_old,kc,k3d)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : qvar, qrho, qu, qv, qw, &
                                   qreint, qpres, version_2, &
                                   npassive, qpass_map, ppm_type, ppm_reference, &
                                   ppm_flatten_before_integrals, &
                                   small_dens, small_pres, gamma_minus_1
    use amrex_constants_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer kc,k3d

    real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt) flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt)   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    real(rt)   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    real(rt) Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)
    real(rt) Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)

    real(rt) qxm (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qxp (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qym (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qyp (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)

    real(rt) dt, a_old

    ! Local variables
    integer i, j
    integer n, ipassive

    real(rt) cc, csq, rho, u, v, w, p, rhoe
    real(rt) rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref

    real(rt) drho, du, dv, dw, dp, drhoe
    real(rt) dup, dvp, dpp
    real(rt) dum, dvm, dpm

    real(rt) enth, alpham, alphap, alpha0r, alpha0e
    real(rt) apright, amright, azrright, azeright
    real(rt) apleft, amleft, azrleft, azeleft
    real(rt) xi, xi1
    real(rt) halfdt

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2
    integer, parameter :: igz = 3

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_3d.f90 :: tracexy_ppm")
    end if

    halfdt = half * dt

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).

    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    ! Version 2 includes the source terms in the jump in 
    ! the velocities that goes through the characteristic projection.

    ! *********************************************************************************************
    ! x-direction
    ! *********************************************************************************************

    ! Trace to left and right edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho  = q(i,j,k3d,qrho)
          u    = q(i,j,k3d,qu)
          v    = q(i,j,k3d,qv)
          w    = q(i,j,k3d,qw)
          p    = q(i,j,k3d,qpres)
          rhoe = q(i,j,k3d,qreint)

          cc = c(i,j,k3d)
          csq = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! ******************************************************************************

          if (i .ge. ilo1) then

             ! Plus state on face i

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                (ppm_reference == 1 .and. u - cc >= 0.0d0) ) then
                 rho_ref = rho
                   u_ref = u
                   v_ref = v
                   w_ref = w
                   p_ref = p
                rhoe_ref = rhoe
             else
                 ! This will be the fastest moving state to the left
                 rho_ref = im(i,j,kc,1,1,qrho)
                   u_ref = im(i,j,kc,1,1,qu)
                   v_ref = im(i,j,kc,1,1,qv)
                   w_ref = im(i,j,kc,1,1,qw)
                   p_ref = im(i,j,kc,1,1,qpres)
                rhoe_ref = im(i,j,kc,1,1,qreint)
             endif
   
             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)
   
             dum   =   u_ref - im(i,j,kc,1,1,qu)
             dpm   =   p_ref - im(i,j,kc,1,1,qpres)
   
             drho  =  rho_ref - im(i,j,kc,1,2,qrho)
             dv    =    v_ref - im(i,j,kc,1,2,qv)
             dw    =    w_ref - im(i,j,kc,1,2,qw)
             dp    =    p_ref - im(i,j,kc,1,2,qpres)
             drhoe = rhoe_ref - im(i,j,kc,1,2,qreint)
   
             dup   =    u_ref - im(i,j,kc,1,3,qu)
             dpp   =    p_ref - im(i,j,kc,1,3,qpres)
   
             if (version_2 .eq. 1) then
                 dum = dum - halfdt*srcq(i,j,k3d,qu)/a_old
                 dup = dup - halfdt*srcq(i,j,k3d,qu)/a_old
             else if (version_2 .eq. 2) then
                 dum = dum - halfdt*im_g(i,j,kc,1,1,igx)/a_old
                 dup = dup - halfdt*im_g(i,j,kc,1,3,igx)/a_old
             end if

            ! These are analogous to the beta's from the original PPM
            ! paper (except we work with rho instead of tau).  This is
            ! simply (l . dq), where dq = qref - I(q)

             alpham = half*(dpm/(rho*cc) - dum)*rho/cc
             alphap = half*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth

             if (u-cc .gt. zero) then
                amright = zero
             else if (u-cc .lt. zero) then
                amright = -alpham
             else
                amright = -half*alpham
             endif

             if (u+cc .gt. zero) then
                apright = zero
             else if (u+cc .lt. zero) then
                apright = -alphap
             else
                apright = -half*alphap
             endif

             if (u .gt. zero) then
                azrright = zero
                azeright = zero
             else if (u .lt. zero) then
                azrright = -alpha0r
                azeright = -alpha0e
             else
                azrright = -half*alpha0r
                azeright = -half*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum(l . dq) r
             qxp(i,j,kc,qrho  ) =  rho_ref +  apright + amright + azrright
             qxp(i,j,kc,qu    ) =    u_ref + (apright - amright)*cc/rho
             qxp(i,j,kc,qreint) = rhoe_ref + (apright + amright)*enth*csq + azeright
             qxp(i,j,kc,qpres ) =    p_ref + (apright + amright)*csq

             ! Transverse velocities -- there's no projection here, so we don't
             ! need a reference state.  We only care about the state traced under
             ! the middle wave
             dv = im(i,j,kc,1,2,qv)
             dw = im(i,j,kc,1,2,qw)
   
             if (version_2 .eq. 1) then
                dv  = dv  - halfdt*srcq(i,j,k3d,qv)/a_old
                dw  = dw  - halfdt*srcq(i,j,k3d,qw)/a_old
             else if (version_2 .eq. 2) then
                dv  = dv  - halfdt*im_g(i,j,kc,1,2,igy)/a_old
                dw  = dw  - halfdt*im_g(i,j,kc,1,2,igz)/a_old
             end if

             ! Recall that I already takes the limit of the parabola
             ! in the event that the wave is not moving toward the
             ! interface
             if (u > zero) then
                qxp(i,j,kc,qv    ) = v
                qxp(i,j,kc,qw    ) = w
             else ! wave moving toward the interface
                qxp(i,j,kc,qv    ) = dv
                qxp(i,j,kc,qw    ) = dw
             endif

             ! We may have already dealt with the flattening in the construction
             ! of the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi  = flatn(i,j,k3d)
                xi1 = one-flatn(i,j,k3d)
 
                qxp(i,j,kc,qrho  ) = xi1*rho  + xi*qxp(i,j,kc,qrho  )
                qxp(i,j,kc,qu    ) = xi1*u    + xi*qxp(i,j,kc,qu    )
                qxp(i,j,kc,qv    ) = xi1*v    + xi*qxp(i,j,kc,qv    )
                qxp(i,j,kc,qw    ) = xi1*w    + xi*qxp(i,j,kc,qw    )
                qxp(i,j,kc,qreint) = xi1*rhoe + xi*qxp(i,j,kc,qreint)
                qxp(i,j,kc,qpres ) = xi1*p    + xi*qxp(i,j,kc,qpres )
             endif

             ! If rho or p too small, set all the slopes to zero
             if (qxp(i,j,kc,qrho ) .lt. small_dens .or. &
                 qxp(i,j,kc,qpres) .lt. small_pres) then
                qxp(i,j,kc,qpres) = p
                qxp(i,j,kc,qrho)  = rho
                qxp(i,j,kc,qu)    = u
             end if

             qxp(i,j,kc,qreint) = qxp(i,j,kc,qpres) / gamma_minus_1

          end if

          ! ******************************************************************************

          if (i .le. ihi1) then

             ! Minus state on face i+1

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                (ppm_reference == 1 .and. u + cc <= 0.0d0) ) then
                 rho_ref = rho
                   u_ref = u
                   v_ref = v
                   w_ref = w
                   p_ref = p
                rhoe_ref = rhoe
             else
                 ! This will be the fastest moving state to the right
                 rho_ref = ip(i,j,kc,1,3,qrho)
                   u_ref = ip(i,j,kc,1,3,qu)
                   v_ref = ip(i,j,kc,1,3,qv)
                   w_ref = ip(i,j,kc,1,3,qw)
                   p_ref = ip(i,j,kc,1,3,qpres)
                rhoe_ref = ip(i,j,kc,1,3,qreint)
             endif
   
             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)

             dum   =    u_ref - ip(i,j,kc,1,1,qu)
             dpm   =    p_ref - ip(i,j,kc,1,1,qpres)
   
             drho  =  rho_ref - ip(i,j,kc,1,2,qrho)
             dp    =    p_ref - ip(i,j,kc,1,2,qpres)
             drhoe = rhoe_ref - ip(i,j,kc,1,2,qreint)

             dup   =    u_ref - ip(i,j,kc,1,3,qu)
             dpp   =    p_ref - ip(i,j,kc,1,3,qpres)
   
             if (version_2 .eq. 1) then
                 dum = dum - halfdt*srcq(i,j,k3d,qu)/a_old
                 dup = dup - halfdt*srcq(i,j,k3d,qu)/a_old
             else if (version_2 .eq. 2) then
                 dum = dum - halfdt*ip_g(i,j,kc,1,1,igx)/a_old
                 dup = dup - halfdt*ip_g(i,j,kc,1,3,igx)/a_old
             end if

             ! These are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This is
             ! simply (l . dq), where dq = qref - I(q)

             alpham = half*(dpm/(rho*cc) - dum)*rho/cc
             alphap = half*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth

             if (u-cc .gt. zero) then
                amleft = -alpham
             else if (u-cc .lt. zero) then
                amleft = zero
             else
                amleft = -half*alpham
             endif

             if (u+cc .gt. zero) then
                apleft = -alphap
             else if (u+cc .lt. zero) then
                apleft = zero 
             else
                apleft = -half*alphap
             endif
   
             if (u .gt. zero) then
                azrleft = -alpha0r
                azeleft = -alpha0e
             else if (u .lt. zero) then
                azrleft = zero
                azeleft = zero
             else
                azrleft = -half*alpha0r
                azeleft = -half*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             qxm(i+1,j,kc,qrho  ) =  rho_ref + apleft + amleft + azrleft
             qxm(i+1,j,kc,qu    ) =    u_ref + (apleft - amleft)*cc/rho
             qxm(i+1,j,kc,qreint) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft
             qxm(i+1,j,kc,qpres ) =    p_ref + (apleft + amleft)*csq

             ! Transverse velocities
             dv    = ip(i,j,kc,1,2,qv)
             dw    = ip(i,j,kc,1,2,qw)
   
             if (version_2 .eq. 1) then
                 dv  = dv  - halfdt*srcq(i,j,k3d,qv)/a_old
                 dw  = dw  - halfdt*srcq(i,j,k3d,qw)/a_old
             else if (version_2 .eq. 2) then
                 dv  = dv  - halfdt*ip_g(i,j,kc,1,2,igy)/a_old
                 dw  = dw  - halfdt*ip_g(i,j,kc,1,2,igz)/a_old
             end if

             if (u < zero) then
                qxm(i+1,j,kc,qv    ) = v
                qxm(i+1,j,kc,qw    ) = w
             else
                qxm(i+1,j,kc,qv    ) = dv
                qxm(i+1,j,kc,qw    ) = dw
             endif
 
             ! We may have already dealt with flattening in the parabolas
             if (ppm_flatten_before_integrals == 0) then
                xi  = flatn(i,j,k3d)
                xi1 = one - flatn(i,j,k3d)
 
                qxm(i+1,j,kc,qrho  ) = xi1*rho  + xi*qxm(i+1,j,kc,qrho  )
                qxm(i+1,j,kc,qu    ) = xi1*u    + xi*qxm(i+1,j,kc,qu    )
                qxm(i+1,j,kc,qv    ) = xi1*v    + xi*qxm(i+1,j,kc,qv    )
                qxm(i+1,j,kc,qw    ) = xi1*w    + xi*qxm(i+1,j,kc,qw    )
                qxm(i+1,j,kc,qreint) = xi1*rhoe + xi*qxm(i+1,j,kc,qreint)
                qxm(i+1,j,kc,qpres ) = xi1*p    + xi*qxm(i+1,j,kc,qpres )
             endif
 
             ! If rho or p too small, set all the slopes to zero
             if (qxm(i+1,j,kc,qrho ) .lt. small_dens .or. &
                 qxm(i+1,j,kc,qpres) .lt. small_pres) then
                qxm(i+1,j,kc,qrho)  = rho
                qxm(i+1,j,kc,qpres) = p
                qxm(i+1,j,kc,qu)    = u
             end if

             qxm(i+1,j,kc,qreint) = qxm(i+1,j,kc,qpres) / gamma_minus_1
          end if

       end do
    end do

    ! ******************************************************************************
    ! Passively advected quantities 
    ! ******************************************************************************

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1

          ! Plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,qu)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = one
             endif
 
             ! The flattening here is a little confusing.  If
             ! ppm_flatten_before_integrals = 0, then we are blending
             ! the cell centered state and the edge state here through
             ! the flattening procedure.  Otherwise, we've already
             ! took care of flattening.  What we want to do is:
             !
             ! q_l*  (1-xi)*q_i + xi*q_l
             !
             ! where
             !
             ! q_l = q_ref - Proj{(q_ref - I)}
             !
             ! and Proj{} represents the characteristic projection.
             ! But for these, there is only 1-wave that matters, the u
             ! wave, so no projection is needed.  Since we are not
             ! projecting, the reference state doesn't matter, so we
             ! take it to be q_i, therefore, we reduce to
             !
             ! q_l* = (1-xi)*q_i + xi*[q_i - (q_i - I)]
             !      = q_i + xi*(I - q_i)
 
             if (u .gt. zero) then
                qxp(i,j,kc,n) = q(i,j,k3d,n)
             else if (u .lt. zero) then
                qxp(i,j,kc,n) = q(i,j,k3d,n) + xi*(im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else
                qxp(i,j,kc,n) = q(i,j,k3d,n) + half*xi*(im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif

      enddo

          ! Minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,qu)
 
             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = one
             endif
 
             if (u .gt. zero) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + xi*(ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else if (u .lt. zero) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n)
             else
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + half*xi*(ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

       enddo
    enddo

    ! *********************************************************************************************
    ! y-direction
    ! *********************************************************************************************

    ! Trace to bottom and top edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho = q(i,j,k3d,qrho)
          u = q(i,j,k3d,qu)
          v = q(i,j,k3d,qv)
          w = q(i,j,k3d,qw)
          p = q(i,j,k3d,qpres)
          rhoe = q(i,j,k3d,qreint)

          cc = c(i,j,k3d)
          csq = cc**2
          enth = ( (rhoe+p)/rho )/csq

          if (j .ge. ilo2) then

             ! Plus state on face j

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                (ppm_reference == 1 .and. v - cc >= 0.0d0) ) then
                 rho_ref = rho
                   u_ref = u
                   v_ref = v
                   w_ref = w
                   p_ref = p
                rhoe_ref = rhoe
             else
                 ! This will be the fastest moving state to the left
                 rho_ref = im(i,j,kc,2,1,qrho)
                   u_ref = im(i,j,kc,2,1,qu)
                   v_ref = im(i,j,kc,2,1,qv)
                   w_ref = im(i,j,kc,2,1,qw)
                   p_ref = im(i,j,kc,2,1,qpres)
                rhoe_ref = im(i,j,kc,2,1,qreint)
             endif
   
             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm   =    v_ref - im(i,j,kc,2,1,qv)
             dpm   =    p_ref - im(i,j,kc,2,1,qpres)
   
             drho  =  rho_ref - im(i,j,kc,2,2,qrho)
             dp    =    p_ref - im(i,j,kc,2,2,qpres)
             drhoe = rhoe_ref - im(i,j,kc,2,2,qreint)

             dvp   =    v_ref - im(i,j,kc,2,3,qv)
             dpp   =    p_ref - im(i,j,kc,2,3,qpres)
   
             if (version_2 .eq. 1) then
                 dvm = dvm - halfdt*srcq(i,j,k3d,qv)/a_old
                 dvp = dvp - halfdt*srcq(i,j,k3d,qv)/a_old
             else if (version_2 .eq. 2) then
                 dvm = dvm - halfdt*im_g(i,j,kc,2,1,igy)/a_old
                 dvp = dvp - halfdt*im_g(i,j,kc,2,3,igy)/a_old
             end if

             ! These are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This
             ! is simply (l . dq), where dq = qref - I(q)
 
             alpham = half*(dpm/(rho*cc) - dvm)*rho/cc
             alphap = half*(dpp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
 
             if (v-cc .gt. zero) then
                amright = zero
             else if (v-cc .lt. zero) then
                amright = -alpham
             else
                amright = -half*alpham
             endif
 
             if (v+cc .gt. zero) then
                apright = zero
             else if (v+cc .lt. zero) then
                apright = -alphap
             else
                apright = -half*alphap
             endif
 
             if (v .gt. zero) then
                azrright = zero
                azeright = zero
             else if (v .lt. zero) then
                azrright = -alpha0r
                azeright = -alpha0e
             else
                azrright = -half*alpha0r
                azeright = -half*alpha0e
             endif
 
             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             qyp(i,j,kc,qrho  ) =  rho_ref +  apright + amright + azrright
             qyp(i,j,kc,qv    ) =    v_ref + (apright - amright)*cc/rho
             qyp(i,j,kc,qreint) = rhoe_ref + (apright + amright)*enth*csq + azeright
             qyp(i,j,kc,qpres ) =    p_ref + (apright + amright)*csq

             ! Transverse velocities
             du    = im(i,j,kc,2,2,qu)
             dw    = im(i,j,kc,2,2,qw)
   
             if (version_2 .eq. 1) then
                 du  = du  - halfdt*srcq(i,j,k3d,qu)/a_old
                 dw  = dw  - halfdt*srcq(i,j,k3d,qw)/a_old
             else if (version_2 .eq. 2) then
                 du  = du  - halfdt*im_g(i,j,kc,2,2,igx)/a_old
                 dw  = dw  - halfdt*im_g(i,j,kc,2,2,igz)/a_old
             end if
 
             if (v > zero) then
                qyp(i,j,kc,qu    ) = u
                qyp(i,j,kc,qw    ) = w
             else ! wave moving toward the interface
                qyp(i,j,kc,qu    ) = du
                qyp(i,j,kc,qw    ) = dw
             endif
 
             ! We may have already dealt with flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi  = flatn(i,j,k3d)
                xi1 = one - flatn(i,j,k3d)
 
                qyp(i,j,kc,qrho  ) = xi1*rho  + xi*qyp(i,j,kc,qrho  )
                qyp(i,j,kc,qv    ) = xi1*v    + xi*qyp(i,j,kc,qv    )
                qyp(i,j,kc,qu    ) = xi1*u    + xi*qyp(i,j,kc,qu    )
                qyp(i,j,kc,qw    ) = xi1*w    + xi*qyp(i,j,kc,qw    )
                qyp(i,j,kc,qreint) = xi1*rhoe + xi*qyp(i,j,kc,qreint)
                qyp(i,j,kc,qpres ) = xi1*p    + xi*qyp(i,j,kc,qpres )
             endif

             ! If rho or p too small, set all the slopes to zero
             if (qyp(i,j,kc,qrho ) .lt. small_dens .or. &
                 qyp(i,j,kc,qpres) .lt. small_pres) then
                qyp(i,j,kc,qrho)  = rho
                qyp(i,j,kc,qpres) = p
                qyp(i,j,kc,qv)    = v
             end if

             qyp(i,j,kc,qreint) = qyp(i,j,kc,qpres) / gamma_minus_1
          end if

          ! ******************************************************************************

          if (j .le. ihi2) then

             ! Minus state on face j+1

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                (ppm_reference == 1 .and. v + cc <= 0.0d0) ) then
                 rho_ref = rho
                   u_ref = u
                   v_ref = v
                   w_ref = w
                   p_ref = p
                rhoe_ref = rhoe
             else
                 ! This will be the fastest moving state to the right
                 rho_ref = ip(i,j,kc,2,3,qrho)
                   u_ref = ip(i,j,kc,2,3,qu)
                   v_ref = ip(i,j,kc,2,3,qv)
                   w_ref = ip(i,j,kc,2,3,qw)
                   p_ref = ip(i,j,kc,2,3,qpres)
                rhoe_ref = ip(i,j,kc,2,3,qreint)
             endif
   
             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm   =    v_ref - ip(i,j,kc,2,1,qv)
             dpm   =    p_ref - ip(i,j,kc,2,1,qpres)
   
             drho  =  rho_ref - ip(i,j,kc,2,2,qrho)
             dp    =    p_ref - ip(i,j,kc,2,2,qpres)
             drhoe = rhoe_ref - ip(i,j,kc,2,2,qreint)

             dvp   =    v_ref - ip(i,j,kc,2,3,qv)
             dpp   =    p_ref - ip(i,j,kc,2,3,qpres)

             if (version_2 .eq. 1) then
                 dvm = dvm - halfdt*srcq(i,j,k3d,qv)/a_old
                 dvp = dvp - halfdt*srcq(i,j,k3d,qv)/a_old
             else if (version_2 .eq. 2) then
                 dvm = dvm - halfdt*ip_g(i,j,kc,2,1,igy)/a_old
                 dvp = dvp - halfdt*ip_g(i,j,kc,2,3,igy)/a_old
             end if

             ! These are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
 
             alpham = half*(dpm/(rho*cc) - dvm)*rho/cc
             alphap = half*(dpp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
 
             if (v-cc .gt. zero) then
                amleft = -alpham
             else if (v-cc .lt. zero) then
                amleft = zero
             else
                amleft = -half*alpham
             endif

             if (v+cc .gt. zero) then
                apleft = -alphap
             else if (v+cc .lt. zero) then
                apleft = zero
             else
                apleft = -half*alphap
             endif

             if (v .gt. zero) then
                azrleft = -alpha0r
                azeleft = -alpha0e
             else if (v .lt. zero) then
                azrleft = zero
                azeleft = zero
             else
                azrleft = -half*alpha0r
                azeleft = -half*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             qym(i,j+1,kc,qrho  ) =  rho_ref +  apleft + amleft + azrleft
             qym(i,j+1,kc,qv    ) =    v_ref + (apleft - amleft)*cc/rho
             qym(i,j+1,kc,qreint) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft
             qym(i,j+1,kc,qpres ) =    p_ref + (apleft + amleft)*csq

             ! Transverse velocities
             du    = ip(i,j,kc,2,2,qu)
             dw    = ip(i,j,kc,2,2,qw)

             if (version_2 .eq. 1) then
                 du  = du  - halfdt*srcq(i,j,k3d,qu)/a_old
                 dw  = dw  - halfdt*srcq(i,j,k3d,qw)/a_old
             else if (version_2 .eq. 2) then
                 du  = du  - halfdt*ip_g(i,j,kc,2,2,igx)/a_old
                 dw  = dw  - halfdt*ip_g(i,j,kc,2,2,igz)/a_old
             end if
 
             if (v < zero) then
                qym(i,j+1,kc,qu    ) = u
                qym(i,j+1,kc,qw    ) = w
             else ! wave is moving toward the interface
                qym(i,j+1,kc,qu    ) = du
                qym(i,j+1,kc,qw    ) = dw
             endif

             ! We may have already dealt with flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi  = flatn(i,j,k3d)
                xi1 = one - flatn(i,j,k3d)
 
                qym(i,j+1,kc,qrho  ) = xi1*rho  + xi*qym(i,j+1,kc,qrho  )
                qym(i,j+1,kc,qv    ) = xi1*v    + xi*qym(i,j+1,kc,qv    )
                qym(i,j+1,kc,qu    ) = xi1*u    + xi*qym(i,j+1,kc,qu    )
                qym(i,j+1,kc,qw    ) = xi1*w    + xi*qym(i,j+1,kc,qw    )
                qym(i,j+1,kc,qreint) = xi1*rhoe + xi*qym(i,j+1,kc,qreint)
                qym(i,j+1,kc,qpres ) = xi1*p    + xi*qym(i,j+1,kc,qpres )
             endif

             ! If rho or p too small, set all the slopes to zero
             if (qym(i,j+1,kc,qrho ) .lt. small_dens .or. &
                 qym(i,j+1,kc,qpres) .lt. small_pres) then
                qym(i,j+1,kc,qrho)  = rho
                qym(i,j+1,kc,qpres) = p
                qym(i,j+1,kc,qv)    = v
             end if

             qym(i,j+1,kc,qreint) = qym(i,j+1,kc,qpres) / gamma_minus_1
          end if

       end do
    end do

    !--------------------------------------------------------------------------
    ! Passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do i = ilo1-1, ihi1+1

          ! Plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,k3d,qv)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = one
             endif

             if (v .gt. zero) then
                qyp(i,j,kc,n) = q(i,j,k3d,n)
             else if (v .lt. zero) then
                qyp(i,j,kc,n) = q(i,j,k3d,n) + xi*(im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else
                qyp(i,j,kc,n) = q(i,j,k3d,n) + half*xi*(im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
          
          ! Minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,qv)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = one
             endif

             if (v .gt. zero) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + xi*(ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else if (v .lt. zero) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n)
             else
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + half*xi*(ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
          
       enddo
    enddo

    end subroutine tracexy_ppm

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

    subroutine tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                          Ip,Im,Ip_g,Im_g, &
                          qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                          srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                          ilo1,ilo2,ihi1,ihi2,dt,a_old,km,kc,k3d)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : qvar, qrho, qu, qv, qw, &
                                   qreint, qpres, version_2, &
                                   npassive, qpass_map, ppm_type, &
                                   npassive, qpass_map, ppm_type, ppm_reference, &
                                   ppm_flatten_before_integrals, &
                                   small_dens, small_pres, gamma_minus_1
    use amrex_constants_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer km,kc,k3d

    real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt) flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt)   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    real(rt)   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    real(rt) Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)
    real(rt) Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)

    real(rt) qzm (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qzp (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    real(rt) dt, a_old

    !     Local variables
    integer i, j
    integer n, ipassive

    real(rt) cc, csq
    real(rt) rho, u, v, w, p, rhoe
    real(rt) rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref
    real(rt) dwp, dpp
    real(rt) dwm, dpm

    real(rt) drho, du, dv, dp, drhoe
    real(rt) enth, alpham, alphap, alpha0r, alpha0e
    real(rt) apright, amright, azrright, azeright
    real(rt) apleft, amleft, azrleft, azeleft
    real(rt) halfdt
    real(rt) xi, xi1

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2
    integer, parameter :: igz = 3

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_3d.f90 :: tracez_ppm")
    end if

    halfdt = half * dt

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! Trace to left and right edges using upwind PPM

    ! Note: in contrast to the above code for x and y, here the loop
    ! is over interfaces, not over cell-centers.

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          ! **************************************************************************
          ! This is all for qzp
          ! **************************************************************************

          rho  = q(i,j,k3d,qrho)
          u    = q(i,j,k3d,qu)
          v    = q(i,j,k3d,qv)
          w    = q(i,j,k3d,qw)
          p    = q(i,j,k3d,qpres)
          rhoe = q(i,j,k3d,qreint)

          cc   = c(i,j,k3d)
          csq  = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! Plus state on face kc

          ! Set the reference state
          if (ppm_reference == 0 .or. &
             (ppm_reference == 1 .and. w - cc >= 0.0d0) ) then
              rho_ref = rho
                u_ref = u
                v_ref = v
                w_ref = w
                p_ref = p
             rhoe_ref = rhoe
          else
                 ! This will be the fastest moving state to the left
              rho_ref = im(i,j,kc,3,1,qrho)
                u_ref = im(i,j,kc,3,1,qu)
                v_ref = im(i,j,kc,3,1,qv)
                w_ref = im(i,j,kc,3,1,qw)
                p_ref = im(i,j,kc,3,1,qpres)
             rhoe_ref = im(i,j,kc,3,1,qreint)
          endif

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   =    w_ref - im(i,j,kc,3,1,qw)
          dpm   =    p_ref - im(i,j,kc,3,1,qpres)

          drho  =  rho_ref - im(i,j,kc,3,2,qrho)
          dp    =    p_ref - im(i,j,kc,3,2,qpres)
          drhoe = rhoe_ref - im(i,j,kc,3,2,qreint)

          dwp   =    w_ref - im(i,j,kc,3,3,qw)
          dpp   =    p_ref - im(i,j,kc,3,3,qpres)

          if (version_2 .eq. 1) then
              dwm = dwm - halfdt*srcq(i,j,k3d,qw)/a_old
              dwp = dwp - halfdt*srcq(i,j,k3d,qw)/a_old
          else if (version_2 .eq. 2) then
              dwm = dwm - halfdt*im_g(i,j,kc,3,1,igz)/a_old
              dwp = dwp - halfdt*im_g(i,j,kc,3,3,igz)/a_old
          end if

          ! These are analogous to the beta's from the original PPM
          ! paper.  This is simply (l . dq), where dq = qref - I(q)
          alpham = half*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = half*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth

          if (w-cc .gt. zero) then
             amright = zero
          else if (w-cc .lt. zero) then
             amright = -alpham
          else
             amright = -half*alpham
          endif
          if (w+cc .gt. zero) then
             apright = zero
          else if (w+cc .lt. zero) then
             apright = -alphap
          else
             apright = -half*alphap
          endif
          if (w .gt. zero) then
             azrright = zero
             azeright = zero
          else if (w .lt. zero) then
             azrright = -alpha0r
             azeright = -alpha0e
          else
             azrright = -half*alpha0r
             azeright = -half*alpha0e
          endif

          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          qzp(i,j,kc,qrho  ) =  rho_ref +  apright + amright + azrright
          qzp(i,j,kc,qw    ) =    w_ref + (apright - amright)*cc/rho
          qzp(i,j,kc,qreint) = rhoe_ref + (apright + amright)*enth*csq + azeright
          qzp(i,j,kc,qpres ) =    p_ref + (apright + amright)*csq

          ! Transverse velocities
          du    = im(i,j,kc,3,2,qu)
          dv    = im(i,j,kc,3,2,qv)

          if (version_2 .eq. 1) then
              du  = du  - halfdt*srcq(i,j,k3d,qu)/a_old
              dv  = dv  - halfdt*srcq(i,j,k3d,qv)/a_old
          else if (version_2 .eq. 2) then
              du  = du  - halfdt*im_g(i,j,kc,3,2,igx)/a_old
              dv  = dv  - halfdt*im_g(i,j,kc,3,2,igy)/a_old
          end if

          if (w > zero) then 
             qzp(i,j,kc,qu    ) = u
             qzp(i,j,kc,qv    ) = v
          else ! wave moving toward the interface
             qzp(i,j,kc,qu    ) = du
             qzp(i,j,kc,qv    ) = dv
          endif

          ! We may have already dealt with flattening in the parabola
          if (ppm_flatten_before_integrals == 0) then
             xi  = flatn(i,j,k3d)
             xi1 = one - flatn(i,j,k3d)

             qzp(i,j,kc,qrho  ) = xi1*rho  + xi*qzp(i,j,kc,qrho  )
             qzp(i,j,kc,qw    ) = xi1*w    + xi*qzp(i,j,kc,qw    )
             qzp(i,j,kc,qu    ) = xi1*u    + xi*qzp(i,j,kc,qu    )
             qzp(i,j,kc,qv    ) = xi1*v    + xi*qzp(i,j,kc,qv    )
             qzp(i,j,kc,qreint) = xi1*rhoe + xi*qzp(i,j,kc,qreint)
             qzp(i,j,kc,qpres ) = xi1*p    + xi*qzp(i,j,kc,qpres )
          endif

          ! If rho or p too small, set all the slopes to zero
         if (qzp(i,j,kc,qrho ) .lt. small_dens .or. &
             qzp(i,j,kc,qpres) .lt. small_pres) then
             qzp(i,j,kc,qrho)  = rho
             qzp(i,j,kc,qpres) = p
             qzp(i,j,kc,qw)    = w
          end if

          qzp(i,j,kc,qreint) = qzp(i,j,kc,qpres) / gamma_minus_1

          ! **************************************************************************
          ! This is all for qzm
          ! **************************************************************************

          ! Minus state on face kc

          ! Note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1

          rho  = q(i,j,k3d-1,qrho)
          u    = q(i,j,k3d-1,qu)
          v    = q(i,j,k3d-1,qv)
          w    = q(i,j,k3d-1,qw)
          p    = q(i,j,k3d-1,qpres)
          rhoe = q(i,j,k3d-1,qreint)

          cc   = c(i,j,k3d-1)
          csq  = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! Set the reference state
          if (ppm_reference == 0 .or. &
             (ppm_reference == 1 .and. w + cc <= 0.0d0) ) then
              rho_ref = rho
                u_ref = u
                v_ref = v
                w_ref = w
                p_ref = p
             rhoe_ref = rhoe
          else
              ! This will be the fastest moving state to the right
              rho_ref = ip(i,j,km,3,3,qrho)
                u_ref = ip(i,j,km,3,3,qu)
                v_ref = ip(i,j,km,3,3,qv)
                w_ref = ip(i,j,km,3,3,qw)
                p_ref = ip(i,j,km,3,3,qpres)
             rhoe_ref = ip(i,j,km,3,3,qreint)
          endif

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   = (   w_ref - ip(i,j,km,3,1,qw))
          dpm   = (   p_ref - ip(i,j,km,3,1,qpres))

          drho  = ( rho_ref - ip(i,j,km,3,2,qrho))
          dp    = (   p_ref - ip(i,j,km,3,2,qpres))
          drhoe = (rhoe_ref - ip(i,j,km,3,2,qreint))

          dwp   = (   w_ref - ip(i,j,km,3,3,qw))
          dpp   = (   p_ref - ip(i,j,km,3,3,qpres))

          if (version_2 .eq. 1) then
              dwm = dwm - halfdt*srcq(i,j,k3d-1,qw)/a_old
              dwp = dwp - halfdt*srcq(i,j,k3d-1,qw)/a_old
          else if (version_2 .eq. 2) then
              dwm = dwm - halfdt*ip_g(i,j,km,3,1,igz)/a_old
              dwp = dwp - halfdt*ip_g(i,j,km,3,3,igz)/a_old
          end if

          ! These are analogous to the beta's from the original PPM
          ! paper.  This is simply (l . dq), where dq = qref - I(q)

          alpham = half*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = half*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
             
          if (w-cc .gt. zero) then
             amleft = -alpham
          else if (w-cc .lt. zero) then
             amleft = zero
          else
             amleft = -half*alpham
          endif
          if (w+cc .gt. zero) then
             apleft = -alphap
          else if (w+cc .lt. zero) then
             apleft = zero
          else
             apleft = -half*alphap
          endif
          if (w .gt. zero) then
             azrleft = -alpha0r
             azeleft = -alpha0e
          else if (w .lt. zero) then
             azrleft = zero
             azeleft = zero
          else
             azrleft = -half*alpha0r
             azeleft = -half*alpha0e
          endif
          
          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          qzm(i,j,kc,qrho  ) =  rho_ref +  apleft + amleft + azrleft
          qzm(i,j,kc,qw    ) =    w_ref + (apleft - amleft)*cc/rho
          qzm(i,j,kc,qreint) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft
          qzm(i,j,kc,qpres ) =    p_ref + (apleft + amleft)*csq

          ! Transverse velocity
          du = ip(i,j,km,3,2,qu)
          dv = ip(i,j,km,3,2,qv)

          if (version_2 .eq. 1) then
              du  = du  - halfdt*srcq(i,j,k3d-1,qu)/a_old
              dv  = dv  - halfdt*srcq(i,j,k3d-1,qv)/a_old
          else if (version_2 .eq. 2) then
              du  = du  - halfdt*ip_g(i,j,km,3,2,igx)/a_old
              dv  = dv  - halfdt*ip_g(i,j,km,3,2,igy)/a_old
          end if
 
          if (w < zero) then
             qzm(i,j,kc,qu    ) = u
             qzm(i,j,kc,qv    ) = v
          else ! wave moving toward the interface
             qzm(i,j,kc,qu    ) = du
             qzm(i,j,kc,qv    ) = dv
          endif

          ! We may have already taken care of flattening in the parabola
          if (ppm_flatten_before_integrals == 0) then
             xi  = flatn(i,j,k3d-1)
             xi1 = one - flatn(i,j,k3d-1)
 
             qzm(i,j,kc,qrho  ) = xi1*rho  + xi*qzm(i,j,kc,qrho  )
             qzm(i,j,kc,qw    ) = xi1*w    + xi*qzm(i,j,kc,qw    )
             qzm(i,j,kc,qu    ) = xi1*u    + xi*qzm(i,j,kc,qu    )
             qzm(i,j,kc,qv    ) = xi1*v    + xi*qzm(i,j,kc,qv    )
             qzm(i,j,kc,qreint) = xi1*rhoe + xi*qzm(i,j,kc,qreint)
             qzm(i,j,kc,qpres ) = xi1*p    + xi*qzm(i,j,kc,qpres )
 
          endif

          ! If rho or p too small, set all the slopes to zero
          if (qzm(i,j,kc,qrho ) .lt. small_dens .or. &
              qzm(i,j,kc,qpres) .lt. small_pres) then
             qzm(i,j,kc,qrho)  = rho
             qzm(i,j,kc,qpres) = p
             qzm(i,j,kc,qw)    = w
          end if

          qzm(i,j,kc,qreint) = qzm(i,j,kc,qpres) / gamma_minus_1

       end do
    end do

    !--------------------------------------------------------------------------
    ! Passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
         n = qpass_map(ipassive)
         do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1

                  ! Plus state on face kc
                  w = q(i,j,k3d,qw)

                  if (ppm_flatten_before_integrals == 0) then
                     xi = flatn(i,j,k3d)
                  else
                     xi = one
                  endif

                  if (w .gt. zero) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n)
                  else if (w .lt. zero) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n) + xi*(im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  else
                     qzp(i,j,kc,n) = q(i,j,k3d,n) + half*xi*(im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  endif

                  ! Minus state on face k
                  w = q(i,j,k3d-1,qw)
                  
                  if (ppm_flatten_before_integrals == 0) then
                     xi = flatn(i,j,k3d-1)
                  else
                     xi = one
                  endif

                  if (w .gt. zero) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) + xi*(ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  else if (w .lt. zero) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n)
                  else
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) + half*xi*(ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  endif

               enddo
         enddo
    enddo

    end subroutine tracez_ppm

end module trace_ppm_module
````

