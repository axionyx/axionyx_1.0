
# File trace\_3d.f90

[**File List**](files.md) **>** [**HydroFortran**](dir_1fab266cd447ad3f3624320661f845f1.md) **>** [**trace\_3d.f90**](trace__3d_8f90.md)

[Go to the documentation of this file.](trace__3d_8f90.md) 


````cpp
! :::
! ::: ------------------------------------------------------------------
! :::
      subroutine tracexy(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dqx,dqy,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
                         qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d,a_old)

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : qvar, qrho, qu, qv, qw, &
                                     qreint, qpres, &
                                     ppm_type, small_dens, small_pres, &
                                     npassive, qpass_map, gamma_minus_1
      use amrex_constants_module
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer kc,k3d

      real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      real(rt)  dqx(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      real(rt)  dqy(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)

      real(rt) qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) a_old
      real(rt) dx, dy, dt

      ! Local variables
      integer i, j, n

      real(rt) dtdx, dtdy
      real(rt) cc, csq, rho, u, v, w, p, rhoe
      real(rt) drho, du, dv, dw, dp, drhoe

      real(rt) enth, alpham, alphap, alpha0r, alpha0e
      real(rt) alpha0u, alpha0v, alpha0w
      real(rt) spminus, spplus, spzero
      real(rt) apright, amright, azrright, azeright
      real(rt) azu1rght, azv1rght, azw1rght
      real(rt) apleft, amleft, azrleft, azeleft
      real(rt) :: azu1left, azv1left, azw1left
      integer          :: ipassive

      dtdx = dt/(dx*a_old)
      dtdy = dt/(dy*a_old)

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in tracexy with ppm_type != 0'
        call amrex_error("Error:: Nyx_advection_3d.f90 :: tracexy")
      end if

      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!

      ! Compute left and right traced states
      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,qrho)
            u = q(i,j,k3d,qu)
            v = q(i,j,k3d,qv)
            w = q(i,j,k3d,qw)
            p = q(i,j,k3d,qpres)
            rhoe = q(i,j,k3d,qreint)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqx(i,j,kc,qrho)
            du = dqx(i,j,kc,qu)
            dv = dqx(i,j,kc,qv)
            dw = dqx(i,j,kc,qw)
            dp = dqx(i,j,kc,qpres)
            drhoe = dqx(i,j,kc,qreint)

            alpham = half*(dp/(rho*cc) - du)*rho/cc
            alphap = half*(dp/(rho*cc) + du)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0v = dv
            alpha0w = dw

            if (u-cc .gt. zero) then
               spminus = -one
            else
               spminus = (u-cc)*dtdx
            endif
            if (u+cc .gt. zero) then
               spplus = -one
            else
               spplus = (u+cc)*dtdx
            endif
            if (u .gt. zero) then
               spzero = -one
            else
               spzero = u*dtdx
            endif

            apright = half*(-one - spplus )*alphap
            amright = half*(-one - spminus)*alpham
            azrright = half*(-one - spzero )*alpha0r
            azeright = half*(-one - spzero )*alpha0e
            azv1rght = half*(-one - spzero )*alpha0v
            azw1rght = half*(-one - spzero )*alpha0w

            if (i .ge. ilo1) then
               qxp(i,j,kc,qrho  ) = rho  + apright + amright + azrright
               qxp(i,j,kc,qu    ) = u    + (apright - amright)*cc/rho
               qxp(i,j,kc,qv    ) = v    +  azv1rght
               qxp(i,j,kc,qw    ) = w    +  azw1rght
               qxp(i,j,kc,qpres ) = p    + (apright + amright)*csq
!              qxp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

               ! If rho or p too small, set all the slopes to zero
               if (qxp(i,j,kc,qrho ) .lt. small_dens .or. &
                   qxp(i,j,kc,qpres) .lt. small_pres) then
                  qxp(i,j,kc,qpres) = p
                  qxp(i,j,kc,qrho)  = rho
                  qxp(i,j,kc,qu)    = u
               end if

               qxp(i,j,kc,qreint) = qxp(i,j,kc,qpres) / gamma_minus_1

            end if

            if (u-cc .ge. zero) then
               spminus = (u-cc)*dtdx
            else
               spminus = one
            endif
            if (u+cc .ge. zero) then
               spplus = (u+cc)*dtdx
            else
               spplus = one
            endif
            if (u .ge. zero) then
               spzero = u*dtdx
            else
               spzero = one
            endif

            apleft = half*(one - spplus )*alphap
            amleft = half*(one - spminus)*alpham
            azrleft = half*(one - spzero )*alpha0r
            azeleft = half*(one - spzero )*alpha0e
            azv1left = half*(one - spzero )*alpha0v
            azw1left = half*(one - spzero )*alpha0w

            if (i .le. ihi1) then
               qxm(i+1,j,kc,qrho  ) = rho  +  apleft + amleft + azrleft
               qxm(i+1,j,kc,qu    ) = u    + (apleft - amleft)*cc/rho
               qxm(i+1,j,kc,qv    ) = v    +  azv1left
               qxm(i+1,j,kc,qw    ) = w    +  azw1left
               qxm(i+1,j,kc,qpres ) = p    + (apleft + amleft)*csq
            !  qxm(i+1,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

               ! If rho or p too small, set all the slopes to zero
               if (qxm(i+1,j,kc,qrho ) .lt. small_dens .or. &
                   qxm(i+1,j,kc,qpres) .lt. small_pres) then
                  qxm(i+1,j,kc,qrho)  = rho
                  qxm(i+1,j,kc,qpres) = p
                  qxm(i+1,j,kc,qu)    = u
               end if
 
               qxm(i+1,j,kc,qreint) = qxm(i+1,j,kc,qpres) / gamma_minus_1

            endif

         enddo
      enddo

      ! Do all of the passively advected quantities in one loop
      do ipassive = 1, npassive
         n = qpass_map(ipassive)

         do j = ilo2-1, ihi2+1
            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,qu)
               if (u .gt. zero) then
                  spzero = -one
               else
                  spzero = u*dtdx
               endif
               qxp(i,j,kc,n) = q(i,j,k3d,n) + half*(-one - spzero )*dqx(i,j,kc,n)
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,qu)
               if (u .ge. zero) then
                  spzero = u*dtdx
               else
                  spzero = one
               endif
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + half*(one - spzero )*dqx(i,j,kc,n)
            enddo
         enddo
      enddo

      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,qrho)
            u = q(i,j,k3d,qu)
            v = q(i,j,k3d,qv)
            w = q(i,j,k3d,qw)
            p = q(i,j,k3d,qpres)
            rhoe = q(i,j,k3d,qreint)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqy(i,j,kc,qrho)
            du = dqy(i,j,kc,qu)
            dv = dqy(i,j,kc,qv)
            dw = dqy(i,j,kc,qw)
            dp = dqy(i,j,kc,qpres)
            drhoe = dqy(i,j,kc,qreint)

            alpham = half*(dp/(rho*cc) - dv)*rho/cc
            alphap = half*(dp/(rho*cc) + dv)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0w = dw

            if (v-cc .gt. zero) then
               spminus = -one
            else
               spminus = (v-cc)*dtdy
            endif
            if (v+cc .gt. zero) then
               spplus = -one
            else
               spplus = (v+cc)*dtdy
            endif
            if (v .gt. zero) then
               spzero = -one
            else
               spzero = v*dtdy
            endif

            apright = half*(-one - spplus )*alphap
            amright = half*(-one - spminus)*alpham
            azrright = half*(-one - spzero )*alpha0r
            azeright = half*(-one - spzero )*alpha0e
            azu1rght = half*(-one - spzero )*alpha0u
            azw1rght = half*(-one - spzero )*alpha0w

            if (j .ge. ilo2) then
               qyp(i,j,kc,qrho) = rho + apright + amright + azrright
               qyp(i,j,kc,qv) = v + (apright - amright)*cc/rho
               qyp(i,j,kc,qu) = u + azu1rght
               qyp(i,j,kc,qw) = w + azw1rght
               qyp(i,j,kc,qpres) = p + (apright + amright)*csq
!              qyp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

               ! If rho or p too small, set all the slopes to zero
               if (qyp(i,j,kc,qrho ) .lt. small_dens .or. &
                   qyp(i,j,kc,qpres) .lt. small_pres) then
                  qyp(i,j,kc,qrho)  = rho
                  qyp(i,j,kc,qpres) = p
                  qyp(i,j,kc,qv)    = v
               end if

               qyp(i,j,kc,qreint) = qyp(i,j,kc,qpres) / gamma_minus_1

            end if

            if (v-cc .ge. zero) then
               spminus = (v-cc)*dtdy
            else
               spminus = one
            endif
            if (v+cc .ge. zero) then
               spplus = (v+cc)*dtdy
            else
               spplus = one
            endif
            if (v .ge. zero) then
               spzero = v*dtdy
            else
               spzero = one
            endif

            apleft = half*(one - spplus )*alphap
            amleft = half*(one - spminus)*alpham
            azrleft = half*(one - spzero )*alpha0r
            azeleft = half*(one - spzero )*alpha0e
            azu1left = half*(one - spzero )*alpha0u
            azw1left = half*(one - spzero )*alpha0w

            if (j .le. ihi2) then
               qym(i,j+1,kc,qrho) = rho + apleft + amleft + azrleft
               qym(i,j+1,kc,qv) = v + (apleft - amleft)*cc/rho
               qym(i,j+1,kc,qu) = u + azu1left
               qym(i,j+1,kc,qw) = w + azw1left
               qym(i,j+1,kc,qpres) = p + (apleft + amleft)*csq
!              qym(i,j+1,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

               ! If rho or p too small, set all the slopes to zero
               if (qym(i,j+1,kc,qrho ) .lt. small_dens .or. &
                   qym(i,j+1,kc,qpres) .lt. small_pres) then
                  qym(i,j+1,kc,qrho)  = rho
                  qym(i,j+1,kc,qpres) = p
                  qym(i,j+1,kc,qv)    = v
               end if

               qym(i,j+1,kc,qreint) = qym(i,j+1,kc,qpres) / gamma_minus_1

            endif

         enddo
      enddo

      ! Do all of the passively advected quantities in one loop
      do ipassive = 1, npassive
         n = qpass_map(ipassive)

         ! Top state
         do j = ilo2, ihi2+1
            do i = ilo1-1, ihi1+1
               v = q(i,j,k3d,qv)
               if (v .gt. zero) then
                  spzero = -one
               else
                  spzero = v*dtdy
               endif
               qyp(i,j,kc,n) = q(i,j,k3d,n) + half*(-one - spzero )*dqy(i,j,kc,n)
            enddo
         end do

         ! Bottom state
         do j = ilo2-1, ihi2
            do i = ilo1-1, ihi1+1
               v = q(i,j,k3d,qv)
               if (v .ge. zero) then
                  spzero = v*dtdy
               else
                  spzero = one
               endif
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + half*(one - spzero )*dqy(i,j,kc,n)
            enddo
         enddo
      enddo

    end subroutine tracexy

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine tracez(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
           dqz,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
           qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
           ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d,a_old)

      use amrex_constants_module
      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : qvar, qrho, qu, qv, qw, &
                                     qreint, qpres, &
                                     ppm_type, small_dens, small_pres, &
                                     npassive, qpass_map, gamma_minus_1

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer km,kc,k3d

      real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      real(rt)  dqz(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      real(rt) qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) a_old
      real(rt) dz, dt

      ! Local variables
      integer i, j
      integer n

      real(rt) dtdz
      real(rt) cc, csq, rho, u, v, w, p, rhoe

      real(rt) drho, du, dv, dw, dp, drhoe
      real(rt) enth, alpham, alphap, alpha0r, alpha0e
      real(rt) alpha0u, alpha0v
      real(rt) spminus, spplus, spzero
      real(rt) apright, amright, azrright, azeright
      real(rt) azu1rght, azv1rght
      real(rt) apleft, amleft, azrleft, azeleft
      real(rt) azu1left, azv1left

      integer ipassive

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in tracez with ppm_type != 0'
        call amrex_error("Error:: Nyx_advection_3d.f90 :: tracez")
      end if

      dtdz = dt/(dz*a_old)

      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!

      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,qrho)
            u = q(i,j,k3d,qu)
            v = q(i,j,k3d,qv)
            w = q(i,j,k3d,qw)
            p = q(i,j,k3d,qpres)
            rhoe = q(i,j,k3d,qreint)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqz(i,j,kc,qrho)
            du = dqz(i,j,kc,qu)
            dv = dqz(i,j,kc,qv)
            dw = dqz(i,j,kc,qw)
            dp = dqz(i,j,kc,qpres)
            drhoe = dqz(i,j,kc,qreint)

            alpham = half*(dp/(rho*cc) - dw)*rho/cc
            alphap = half*(dp/(rho*cc) + dw)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .gt. zero) then
               spminus = -one
            else
               spminus = (w-cc)*dtdz
            endif
            if (w+cc .gt. zero) then
               spplus = -one
            else
               spplus = (w+cc)*dtdz
            endif
            if (w .gt. zero) then
               spzero = -one
            else
               spzero = w*dtdz
            endif

            apright = half*(-one - spplus )*alphap
            amright = half*(-one - spminus)*alpham
            azrright = half*(-one - spzero )*alpha0r
            azeright = half*(-one - spzero )*alpha0e
            azu1rght = half*(-one - spzero )*alpha0u
            azv1rght = half*(-one - spzero )*alpha0v

            qzp(i,j,kc,qrho) = rho + apright + amright + azrright
            qzp(i,j,kc,qw) = w + (apright - amright)*cc/rho
            qzp(i,j,kc,qu) = u + azu1rght
            qzp(i,j,kc,qv) = v + azv1rght
            qzp(i,j,kc,qpres) = p + (apright + amright)*csq
!           qzp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

            ! If rho or p too small, set all the slopes to zero
            if (qzp(i,j,kc,qrho ) .lt. small_dens .or. &
                qzp(i,j,kc,qpres) .lt. small_pres) then
               qzp(i,j,kc,qrho)  = rho
               qzp(i,j,kc,qpres) = p
               qzp(i,j,kc,qw)    = w
            end if
 
            qzp(i,j,kc,qreint) = qzp(i,j,kc,qpres) / gamma_minus_1

            ! **************************************************************************

            ! repeat above with km (k3d-1) to get qzm at kc
            cc  = c(i,j,k3d-1)
            csq = cc**2
            rho = q(i,j,k3d-1,qrho)
            u = q(i,j,k3d-1,qu)
            v = q(i,j,k3d-1,qv)
            w = q(i,j,k3d-1,qw)
            p = q(i,j,k3d-1,qpres)
            rhoe = q(i,j,k3d-1,qreint)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqz(i,j,km,qrho)
            du = dqz(i,j,km,qu)
            dv = dqz(i,j,km,qv)
            dw = dqz(i,j,km,qw)
            dp = dqz(i,j,km,qpres)
            drhoe = dqz(i,j,km,qreint)

            alpham = half*(dp/(rho*cc) - dw)*rho/cc
            alphap = half*(dp/(rho*cc) + dw)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .ge. zero) then
               spminus = (w-cc)*dtdz
            else
               spminus = one
            endif
            if (w+cc .ge. zero) then
               spplus = (w+cc)*dtdz
            else
               spplus = one
            endif
            if (w .ge. zero) then
               spzero = w*dtdz
            else
               spzero = one
            endif

            apleft = half*(one - spplus )*alphap
            amleft = half*(one - spminus)*alpham
            azrleft = half*(one - spzero )*alpha0r
            azeleft = half*(one - spzero )*alpha0e
            azu1left = half*(one - spzero )*alpha0u
            azv1left = half*(one - spzero )*alpha0v

            qzm(i,j,kc,qrho) = rho + apleft + amleft + azrleft
            qzm(i,j,kc,qw) = w + (apleft - amleft)*cc/rho
            qzm(i,j,kc,qu) = u + azu1left
            qzm(i,j,kc,qv) = v + azv1left
            qzm(i,j,kc,qpres) = p + (apleft + amleft)*csq
!           qzm(i,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

            ! If rho or p too small, set all the slopes to zero
            if (qzm(i,j,kc,qrho ) .lt. small_dens .or. &
                qzm(i,j,kc,qpres) .lt. small_pres) then
               qzm(i,j,kc,qrho)  = rho
               qzm(i,j,kc,qpres) = p
               qzm(i,j,kc,qw)    = w
            end if
 
            qzm(i,j,kc,qreint) = qzm(i,j,kc,qpres) / gamma_minus_1

         enddo
      enddo

      ! Do all of the passively advected quantities in one loop
      do ipassive = 1, npassive
         n = qpass_map(ipassive)

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               ! Top state
               w = q(i,j,k3d,qw)
               if (w .gt. zero) then
                  spzero = -one
               else
                  spzero = w*dtdz
               endif
               qzp(i,j,kc,n) = q(i,j,k3d,n) + half*(-one - spzero )*dqz(i,j,kc,n)

               ! Bottom state
               w = q(i,j,k3d-1,qw)
               if (w .ge. zero) then
                  spzero = w*dtdz
               else
                  spzero = one
               endif
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + half*(one - spzero )*dqz(i,j,km,n)

            enddo
         enddo
      enddo

    end subroutine tracez
````

