
# File trans\_3d.f90

[**File List**](files.md) **>** [**HydroFortran**](dir_1fab266cd447ad3f3624320661f845f1.md) **>** [**trans\_3d.f90**](trans__3d_8f90.md)

[Go to the documentation of this file.](trans__3d_8f90.md) 


````cpp
module transverse_module
 
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
 
  implicit none
 
contains

      subroutine transx1(qym,qymo,qyp,qypo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fx,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         cdtdx,ilo,ihi,jlo,jhi,kc,k3d)

      ! Note that what we call ilo here is ilo = lo(1)
      ! Note that what we call ihi here is ihi = hi(1)
      ! Note that what we call jlo here is jlo = lo(2) - 1
      ! Note that what we call jhi here is jhi = hi(2) + 1

      use meth_params_module, only : qvar, nvar, qrho, qu, qv, qw, &
                                     qpres, qreint, &
                                     urho, umx, umy, umz, ueden, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer ilo,ihi,jlo,jhi,kc,k3d

      real(rt)  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      real(rt) ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      real(rt) pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)

      ! Note that cdtdx = dtdx/3.d0/a_half
      real(rt) cdtdx

      integer i, j
      integer n, nq

      real(rt) rrnew, rr
      real(rt) rrry, rrly
      real(rt) rury, ruly
      real(rt) rvry, rvly
      real(rt) rwry, rwly
      real(rt) ekenry, ekenly
      real(rt) pnewry, pnewly
      real(rt) rery, rely
      real(rt) rrnewry, rrnewly
      real(rt) runewry, runewly
      real(rt) rvnewry, rvnewly
      real(rt) rwnewry, rwnewly
      real(rt) renewry, renewly
      real(rt) rhoekenry, rhoekenly
      real(rt) compn, compu
      real(rt) pgp, pgm, ugp, ugm, dup, pav, du

      integer          :: ipassive

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            !DIR$ vector always
            do i = ilo, ihi

               compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n)) 

               if (j.ge.jlo+1) then
                  rr = qyp(i,j,kc,qrho)
                  rrnew = rr - cdtdx*(fx(i+1,j,kc,urho) - fx(i,j,kc,urho))
                  compu = rr*qyp(i,j,kc,nq) - compn
                  qypo(i,j,kc,nq) = compu/rrnew
               end if

               if (j.le.jhi-1) then
                  rr = qym(i,j+1,kc,qrho)
                  rrnew = rr - cdtdx*(fx(i+1,j,kc,urho) - fx(i,j,kc,urho))
                  compu = rr*qym(i,j+1,kc,nq) - compn
                  qymo(i,j+1,kc,nq) = compu/rrnew
               end if

            enddo
         enddo
      enddo

      ! NOTE: it is better *not* to protect against small density in this routine

      do j = jlo, jhi
         do i = ilo, ihi

            pgp = pgdnvx(i+1,j,kc)
            pgm = pgdnvx(i,j,kc)
            ugp = ugdnvx(i+1,j,kc)
            ugm = ugdnvx(i,j,kc)

            if (j.ge.jlo+1) then
               ! Convert to conservation form
               rrry = qyp(i,j,kc,qrho)
               rury = rrry*qyp(i,j,kc,qu)
               rvry = rrry*qyp(i,j,kc,qv)
               rwry = rrry*qyp(i,j,kc,qw)
               ekenry = half*rrry*(qyp(i,j,kc,qu)**2 + qyp(i,j,kc,qv)**2 + qyp(i,j,kc,qw)**2)
               rery = qyp(i,j,kc,qreint) + ekenry

               ! Add transverse terms
               rrnewry = rrry - cdtdx*(fx(i+1,j,kc,urho ) - fx(i,j,kc,urho ))
               runewry = rury - cdtdx*(fx(i+1,j,kc,umx  ) - fx(i,j,kc,umx  ))
               rvnewry = rvry - cdtdx*(fx(i+1,j,kc,umy  ) - fx(i,j,kc,umy  ))
               rwnewry = rwry - cdtdx*(fx(i+1,j,kc,umz  ) - fx(i,j,kc,umz  ))
               renewry = rery - cdtdx*(fx(i+1,j,kc,ueden) - fx(i,j,kc,ueden))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewry .lt. zero) then
                 rrnewry = rrry
                 runewry = rury
                 rvnewry = rvry
                 rwnewry = rwry
                 renewry = rery
               end if
            end if

            if (j.le.jhi-1) then
               rrly = qym(i,j+1,kc,qrho)
               ruly = rrly*qym(i,j+1,kc,qu)
               rvly = rrly*qym(i,j+1,kc,qv)
               rwly = rrly*qym(i,j+1,kc,qw)
               ekenly = half*rrly* &
                    (qym(i,j+1,kc,qu)**2 + qym(i,j+1,kc,qv)**2 + qym(i,j+1,kc,qw)**2)
               rely = qym(i,j+1,kc,qreint) + ekenly

               ! Add transverse terms
               rrnewly = rrly - cdtdx*(fx(i+1,j,kc,urho ) - fx(i,j,kc,urho ))
               runewly = ruly - cdtdx*(fx(i+1,j,kc,umx  ) - fx(i,j,kc,umx  ))
               rvnewly = rvly - cdtdx*(fx(i+1,j,kc,umy  ) - fx(i,j,kc,umy  ))
               rwnewly = rwly - cdtdx*(fx(i+1,j,kc,umz  ) - fx(i,j,kc,umz  ))
               renewly = rely - cdtdx*(fx(i+1,j,kc,ueden) - fx(i,j,kc,ueden))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewly .lt. zero) then
                 rrnewly = rrly
                 runewly = ruly
                 rvnewly = rvly
                 rwnewly = rwly
                 renewly = rely
               end if
            end if

            dup = pgp*ugp - pgm*ugm
            pav = half*(pgp+pgm)
            du = ugp-ugm

            ! Convert back to non-conservation form
            if (j.ge.jlo+1) then
               qypo(i,j,kc,qrho) = rrnewry
               qypo(i,j,kc,qu) = runewry/qypo(i,j,kc,qrho)
               qypo(i,j,kc,qv) = rvnewry/qypo(i,j,kc,qrho)
               qypo(i,j,kc,qw) = rwnewry/qypo(i,j,kc,qrho)
               rhoekenry = half*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,kc,qrho)

               qypo(i,j,kc,qreint) = renewry - rhoekenry
               qypo(i,j,kc,qpres) = qypo(i,j,kc,qreint) * gamma_minus_1

               if (qypo(i,j,kc,qpres) .lt. small_pres) then
                   pnewry = qyp(i,j  ,kc,qpres) - cdtdx*(dup + pav*du*gamma_minus_1)
                   qypo(i,j,kc,qpres ) = pnewry
                   qypo(i,j,kc,qreint) = qypo(i,j,kc,qpres) / gamma_minus_1
               end if
            end if

            if (j.le.jhi-1) then
               qymo(i,j+1,kc,qrho) = rrnewly
               qymo(i,j+1,kc,qu) = runewly/qymo(i,j+1,kc,qrho)
               qymo(i,j+1,kc,qv) = rvnewly/qymo(i,j+1,kc,qrho)
               qymo(i,j+1,kc,qw) = rwnewly/qymo(i,j+1,kc,qrho)
               rhoekenly = half*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,kc,qrho)

               qymo(i,j+1,kc,qreint) = renewly - rhoekenly
               qymo(i,j+1,kc,qpres) = qymo(i,j+1,kc,qreint) * gamma_minus_1

               if (qymo(i,j+1,kc,qpres) .lt. small_pres) then
                   pnewly = qym(i,j+1,kc,qpres) - cdtdx*(dup + pav*du*gamma_minus_1)
                   qymo(i,j+1,kc,qpres ) = pnewly
                   qymo(i,j+1,kc,qreint) = qymo(i,j+1,kc,qpres) / gamma_minus_1
               end if
            end if

         enddo
      enddo

      end subroutine transx1

      subroutine transx2(qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fx,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         cdtdx,ilo,ihi,jlo,jhi,kc,km,k3d)

      use meth_params_module, only : qvar, nvar, qrho, qu, qv, qw, &
                                     qpres, qreint, &
                                     urho, umx, umy, umz, ueden, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer ilo,ihi,jlo,jhi,kc,km,k3d

      real(rt)  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      real(rt) ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      real(rt) pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)

      ! Note that cdtdx = dtdx/3.d0/a_half
      real(rt) cdtdx

      integer i, j
      integer n, nq

      real(rt) rrnew, rr
      real(rt) rrrz, rrlz
      real(rt) rurz, rulz
      real(rt) rvrz, rvlz
      real(rt) rwrz, rwlz
      real(rt) ekenrz, ekenlz
      real(rt) rerz, relz
      real(rt) pnewrz, pnewlz
      real(rt) rrnewrz, rrnewlz
      real(rt) runewrz, runewlz
      real(rt) rvnewrz, rvnewlz
      real(rt) rwnewrz, rwnewlz
      real(rt) renewrz, renewlz
      real(rt) rhoekenrz, rhoekenlz
      real(rt) compn, compu
      real(rt) pgp, pgm, ugp, ugm, dup, pav, du

      integer          :: ipassive
      
      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            !DIR$ vector always
            do i = ilo, ihi

                compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

                rr = qzp(i,j,kc,qrho)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,urho) - fx(i,j,kc,urho))
                compu = rr*qzp(i,j,kc,nq) - compn
                qzpo(i,j,kc,nq) = compu/rrnew

                compn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))

                rr = qzm(i,j,kc,qrho)
                rrnew = rr - cdtdx*(fx(i+1,j,km,urho) - fx(i,j,km,urho))
                compu = rr*qzm(i,j,kc,nq) - compn
                qzmo(i,j,kc,nq) = compu/rrnew

            enddo
         enddo
      enddo

      do j = jlo, jhi
          do i = ilo, ihi

             ! ************************************************************************
             ! Convert to conservation form
             rrrz =      qzp(i,j,kc,qrho)
             rurz = rrrz*qzp(i,j,kc,qu)
             rvrz = rrrz*qzp(i,j,kc,qv)
             rwrz = rrrz*qzp(i,j,kc,qw)
             ekenrz = half*rrrz*(qzp(i,j,kc,qu)**2 + qzp(i,j,kc,qv)**2 + qzp(i,j,kc,qw)**2)
             rerz = qzp(i,j,kc,qreint) + ekenrz

             ! Add transverse terms
             rrnewrz = rrrz - cdtdx*(fx(i+1,j,kc,urho ) - fx(i,j,kc,urho ))
             runewrz = rurz - cdtdx*(fx(i+1,j,kc,umx  ) - fx(i,j,kc,umx  ))
             rvnewrz = rvrz - cdtdx*(fx(i+1,j,kc,umy  ) - fx(i,j,kc,umy  ))
             rwnewrz = rwrz - cdtdx*(fx(i+1,j,kc,umz  ) - fx(i,j,kc,umz  ))
             renewrz = rerz - cdtdx*(fx(i+1,j,kc,ueden) - fx(i,j,kc,ueden))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewrz .lt. zero) then
                 rrnewrz = rrrz
                 runewrz = rurz
                 rvnewrz = rvrz
                 rwnewrz = rwrz
                 renewrz = rerz
            end if

             ! Convert back to non-conservation form
             qzpo(i,j,kc,qrho) = rrnewrz
             qzpo(i,j,kc,qu) = runewrz/qzpo(i,j,kc,qrho)
             qzpo(i,j,kc,qv) = rvnewrz/qzpo(i,j,kc,qrho)
             qzpo(i,j,kc,qw) = rwnewrz/qzpo(i,j,kc,qrho)
             rhoekenrz = half*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,qrho)

             qzpo(i,j,kc,qreint) = renewrz - rhoekenrz
             qzpo(i,j,kc,qpres) = qzpo(i,j,kc,qreint) * gamma_minus_1

             if (qzpo(i,j,kc,qpres) .lt. small_pres) then
                 pgp = pgdnvx(i+1,j,kc)
                 pgm = pgdnvx(i,j,kc)
                 ugp = ugdnvx(i+1,j,kc)
                 ugm = ugdnvx(i,j,kc)
                 dup = pgp*ugp - pgm*ugm
                 pav = half*(pgp+pgm)
                 du = ugp-ugm
                 pnewrz = qzp(i,j,kc,qpres) - cdtdx*(dup + pav*du*gamma_minus_1)
                 qzpo(i,j,kc,qpres ) = pnewrz
                 qzpo(i,j,kc,qreint) = qzpo(i,j,kc,qpres) / gamma_minus_1
             end if
             ! ************************************************************************

             ! ************************************************************************
             ! Convert to conservation form
             rrlz =      qzm(i,j,kc,qrho)
             rulz = rrlz*qzm(i,j,kc,qu)
             rvlz = rrlz*qzm(i,j,kc,qv)
             rwlz = rrlz*qzm(i,j,kc,qw)
             ekenlz = half*rrlz*(qzm(i,j,kc,qu)**2 + qzm(i,j,kc,qv)**2 + qzm(i,j,kc,qw)**2)
             relz = qzm(i,j,kc,qreint) + ekenlz

             ! Add transverse terms
             rrnewlz = rrlz - cdtdx*(fx(i+1,j,km,urho ) - fx(i,j,km,urho ))
             runewlz = rulz - cdtdx*(fx(i+1,j,km,umx  ) - fx(i,j,km,umx  ))
             rvnewlz = rvlz - cdtdx*(fx(i+1,j,km,umy  ) - fx(i,j,km,umy  ))
             rwnewlz = rwlz - cdtdx*(fx(i+1,j,km,umz  ) - fx(i,j,km,umz  ))
             renewlz = relz - cdtdx*(fx(i+1,j,km,ueden) - fx(i,j,km,ueden))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewlz .lt. zero) then
                 rrnewlz = rrlz
                 runewlz = rulz
                 rvnewlz = rvlz
                 rwnewlz = rwlz
                 renewlz = relz
            end if

             ! Convert back to non-conservation form
             qzmo(i,j,kc,qrho) = rrnewlz
             qzmo(i,j,kc,qu) = runewlz/qzmo(i,j,kc,qrho)
             qzmo(i,j,kc,qv) = rvnewlz/qzmo(i,j,kc,qrho)
             qzmo(i,j,kc,qw) = rwnewlz/qzmo(i,j,kc,qrho)
             rhoekenlz = half*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,qrho)

             qzmo(i,j,kc,qreint) = renewlz - rhoekenlz
             qzmo(i,j,kc,qpres) = qzmo(i,j,kc,qreint) * gamma_minus_1

             if (qzmo(i,j,kc,qpres) .lt. small_pres) then
                 pgp = pgdnvx(i+1,j,km)
                 pgm = pgdnvx(i,j,km)
                 ugp = ugdnvx(i+1,j,km)
                 ugm = ugdnvx(i,j,km)
                 dup = pgp*ugp - pgm*ugm
                 pav = half*(pgp+pgm)
                 du = ugp-ugm
                 pnewlz = qzm(i,j,kc,qpres) - cdtdx*(dup + pav*du*gamma_minus_1)
                 qzmo(i,j,kc,qpres ) = pnewlz
                 qzmo(i,j,kc,qreint) = qzmo(i,j,kc,qpres) / gamma_minus_1
             end if
             ! ************************************************************************

          enddo
      enddo

      end subroutine transx2

      subroutine transy1(qxm,qxmo,qxp,qxpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fy,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         cdtdy,ilo,ihi,jlo,jhi,kc,k3d)

      use meth_params_module, only : qvar, nvar, qrho, qu, qv, qw, &
                                     qpres, qreint, &
                                     urho, umx, umy, umz, ueden, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer ilo,ihi,jlo,jhi,kc,k3d

      real(rt)  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      real(rt) ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      real(rt) pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)

      ! Note that cdtdx = dtdx/3.d0/a_half
      real(rt) cdtdy

      integer i, j
      integer n, nq

      real(rt) rrnew, rr
      real(rt) compn, compu
      real(rt) rrrx, rrlx
      real(rt) rurx, rulx
      real(rt) rvrx, rvlx
      real(rt) rwrx, rwlx
      real(rt) ekenrx, ekenlx
      real(rt) rerx, relx
      real(rt) rrnewrx, rrnewlx
      real(rt) runewrx, runewlx
      real(rt) rvnewrx, rvnewlx
      real(rt) rwnewrx, rwnewlx
      real(rt) renewrx, renewlx
      real(rt) pnewrx, pnewlx
      real(rt) rhoekenrx, rhoekenlx
      real(rt) pgp, pgm, ugp, ugm
      real(rt) du,dup,pav

      integer          :: ipassive
      
      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            !DIR$ vector always
            do i = ilo, ihi

               compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               if (i.ge.ilo+1) then
                  rr = qxp(i,j,kc,qrho)
                  rrnew = rr - cdtdy*(fy(i,j+1,kc,urho) - fy(i,j,kc,urho))
                  compu = rr*qxp(i,j,kc,nq) - compn
                  qxpo(i,j,kc,nq) = compu/rrnew
               end if

               if (i.le.ihi-1) then
                  rr = qxm(i+1,j,kc,qrho)
                  rrnew = rr - cdtdy*(fy(i,j+1,kc,urho) - fy(i,j,kc,urho))
                  compu = rr*qxm(i+1,j,kc,nq) - compn
                  qxmo(i+1,j,kc,nq) = compu/rrnew
               end if

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

            pgp = pgdnvy(i,j+1,kc)
            pgm = pgdnvy(i,j,kc)
            ugp = ugdnvy(i,j+1,kc)
            ugm = ugdnvy(i,j,kc)

            ! Convert to conservation form
            if (i.ge.ilo+1) then
               rrrx = qxp(i,j,kc,qrho)
               rurx = rrrx*qxp(i,j,kc,qu)
               rvrx = rrrx*qxp(i,j,kc,qv)
               rwrx = rrrx*qxp(i,j,kc,qw)
               ekenrx = half*rrrx*(qxp(i,j,kc,qu)**2 + qxp(i,j,kc,qv)**2 &
                    + qxp(i,j,kc,qw)**2)
               rerx = qxp(i,j,kc,qreint) + ekenrx

               ! Add transverse terms
               rrnewrx = rrrx - cdtdy*(fy(i,j+1,kc,urho ) - fy(i,j,kc,urho ))
               runewrx = rurx - cdtdy*(fy(i,j+1,kc,umx  ) - fy(i,j,kc,umx  ))
               rvnewrx = rvrx - cdtdy*(fy(i,j+1,kc,umy  ) - fy(i,j,kc,umy  ))
               rwnewrx = rwrx - cdtdy*(fy(i,j+1,kc,umz  ) - fy(i,j,kc,umz  ))
               renewrx = rerx - cdtdy*(fy(i,j+1,kc,ueden) - fy(i,j,kc,ueden))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewrx .lt. zero) then
                 rrnewrx = rrrx
                 runewrx = rurx
                 rvnewrx = rvrx
                 rwnewrx = rwrx
                 renewrx = rerx
               end if
            end if
   
            if (i.le.ihi-1) then
               rrlx = qxm(i+1,j,kc,qrho)
               rulx = rrlx*qxm(i+1,j,kc,qu)
               rvlx = rrlx*qxm(i+1,j,kc,qv)
               rwlx = rrlx*qxm(i+1,j,kc,qw)
               ekenlx = half*rrlx*(qxm(i+1,j,kc,qu)**2 + qxm(i+1,j,kc,qv)**2 &
                    + qxm(i+1,j,kc,qw)**2)
               relx = qxm(i+1,j,kc,qreint) + ekenlx

               ! Add transverse terms
               rrnewlx = rrlx - cdtdy*(fy(i,j+1,kc,urho ) - fy(i,j,kc,urho ))
               runewlx = rulx - cdtdy*(fy(i,j+1,kc,umx  ) - fy(i,j,kc,umx  ))
               rvnewlx = rvlx - cdtdy*(fy(i,j+1,kc,umy  ) - fy(i,j,kc,umy  ))
               rwnewlx = rwlx - cdtdy*(fy(i,j+1,kc,umz  ) - fy(i,j,kc,umz  ))
               renewlx = relx - cdtdy*(fy(i,j+1,kc,ueden) - fy(i,j,kc,ueden))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewlx .lt. zero) then
                    rrnewlx = rrlx
                    runewlx = rulx
                    rvnewlx = rvlx
                    rwnewlx = rwlx
                    renewlx = relx
               end if
            end if

            dup = pgp*ugp - pgm*ugm
            pav = half*(pgp+pgm)
            du = ugp-ugm

            ! Convert back to non-conservation form

            ! ************************************************************************
            if (i.ge.ilo+1) then
               qxpo(i,j,kc,qrho) = rrnewrx
               qxpo(i,j,kc,qu) = runewrx/qxpo(i,j,kc,qrho)
               qxpo(i,j,kc,qv) = rvnewrx/qxpo(i,j,kc,qrho)
               qxpo(i,j,kc,qw) = rwnewrx/qxpo(i,j,kc,qrho)
               rhoekenrx = half*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,kc,qrho)

               qxpo(i,j,kc,qreint)= renewrx - rhoekenrx
               qxpo(i,j,kc,qpres) = qxpo(i,j,kc,qreint) * gamma_minus_1

               if (qxpo(i,j,kc,qpres) .lt. small_pres) then
                   pnewrx = qxp(i  ,j,kc,qpres) - cdtdy*(dup + pav*du*gamma_minus_1)
                   qxpo(i,j,kc,qpres) = pnewrx
                   qxpo(i,j,kc,qreint) = qxpo(i,j,kc,qpres) / gamma_minus_1
               end if
            end if
            ! ************************************************************************

            ! ************************************************************************
            if (i.le.ihi-1) then
               qxmo(i+1,j,kc,qrho) = rrnewlx
               qxmo(i+1,j,kc,qu) = runewlx/qxmo(i+1,j,kc,qrho)
               qxmo(i+1,j,kc,qv) = rvnewlx/qxmo(i+1,j,kc,qrho)
               qxmo(i+1,j,kc,qw) = rwnewlx/qxmo(i+1,j,kc,qrho)
               rhoekenlx = half*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,kc,qrho)

               qxmo(i+1,j,kc,qreint)= renewlx - rhoekenlx
               qxmo(i+1,j,kc,qpres) = qxmo(i+1,j,kc,qreint) * gamma_minus_1

               if (qxmo(i+1,j,kc,qpres) .lt. small_pres) then
                   pnewlx = qxm(i+1,j,kc,qpres) - cdtdy*(dup + pav*du*gamma_minus_1)
                   qxmo(i+1,j,kc,qpres ) = pnewlx
                   qxmo(i+1,j,kc,qreint) = qxmo(i+1,j,kc,qpres) / gamma_minus_1
               end if
            end if
            ! ************************************************************************

         enddo
      enddo

      end subroutine transy1

      subroutine transy2(qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fy,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)

      use meth_params_module, only : qvar, nvar, qrho, qu, qv, qw, &
                                     qpres, qreint, &
                                     urho, umx, umy, umz, ueden, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer ilo,ihi,jlo,jhi,kc,km,k3d

      real(rt)  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      real(rt) ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      real(rt) pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)

      ! Note that cdtdy = dtdy/3.d0/a_half
      real(rt) cdtdy

      integer i, j
      integer n, nq

      real(rt) rrnew, rr
      real(rt) compn, compu
      real(rt) rrrz, rrlz
      real(rt) rurz, rulz
      real(rt) rvrz, rvlz
      real(rt) rwrz, rwlz
      real(rt) ekenrz, ekenlz
      real(rt) rerz, relz
      real(rt) rrnewrz, rrnewlz
      real(rt) runewrz, runewlz
      real(rt) rvnewrz, rvnewlz
      real(rt) rwnewrz, rwnewlz
      real(rt) renewrz, renewlz
      real(rt) pnewrz , pnewlz
      real(rt) rhoekenrz, rhoekenlz
      real(rt) pgp, pgm, ugp, ugm, dup, pav, du

      integer          :: ipassive

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            !DIR$ vector always
            do i = ilo, ihi

               compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               rr = qzp(i,j,kc,qrho)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,urho) - fy(i,j,kc,urho))
               compu = rr*qzp(i,j,kc,nq) - compn
               qzpo(i,j,kc,nq) = compu/rrnew

               compn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

               rr = qzm(i,j,kc,qrho)
               rrnew = rr - cdtdy*(fy(i,j+1,km,urho) - fy(i,j,km,urho))
               compu = rr*qzm(i,j,kc,nq) - compn
               qzmo(i,j,kc,nq) = compu/rrnew

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

            ! Convert to conservation form

            ! ************************************************************************
            rrrz = qzp(i,j,kc,qrho)
            rurz = rrrz*qzp(i,j,kc,qu)
            rvrz = rrrz*qzp(i,j,kc,qv)
            rwrz = rrrz*qzp(i,j,kc,qw)
            ekenrz = half*rrrz*(qzp(i,j,kc,qu)**2 + qzp(i,j,kc,qv)**2 &
                 + qzp(i,j,kc,qw)**2)
            rerz = qzp(i,j,kc,qreint) + ekenrz

            ! Add transverse terms
            rrnewrz = rrrz - cdtdy*(fy(i,j+1,kc,urho ) - fy(i,j,kc,urho ))
            runewrz = rurz - cdtdy*(fy(i,j+1,kc,umx  ) - fy(i,j,kc,umx  ))
            rvnewrz = rvrz - cdtdy*(fy(i,j+1,kc,umy  ) - fy(i,j,kc,umy  ))
            rwnewrz = rwrz - cdtdy*(fy(i,j+1,kc,umz  ) - fy(i,j,kc,umz  ))
            renewrz = rerz - cdtdy*(fy(i,j+1,kc,ueden) - fy(i,j,kc,ueden))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewrz .lt. zero) then
                 rrnewrz = rrrz
                 runewrz = rurz
                 rvnewrz = rvrz
                 rwnewrz = rwrz
                 renewrz = rerz
            end if

            ! Convert back to non-conservation form
            qzpo(i,j,kc,qrho) = rrnewrz
            qzpo(i,j,kc,qu) = runewrz/qzpo(i,j,kc,qrho)
            qzpo(i,j,kc,qv) = rvnewrz/qzpo(i,j,kc,qrho)
            qzpo(i,j,kc,qw) = rwnewrz/qzpo(i,j,kc,qrho)
            rhoekenrz = half*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,qrho)

            qzpo(i,j,kc,qreint)= renewrz - rhoekenrz
            qzpo(i,j,kc,qpres) = qzpo(i,j,kc,qreint) * gamma_minus_1

            if (qzpo(i,j,kc,qpres) .lt. small_pres) then
                pgp = pgdnvy(i,j+1,kc)
                pgm = pgdnvy(i,j,kc)
                ugp = ugdnvy(i,j+1,kc)
                ugm = ugdnvy(i,j,kc)
                dup = pgp*ugp - pgm*ugm
                pav = half*(pgp+pgm)
                du = ugp-ugm
                pnewrz = qzp(i,j,kc,qpres) - cdtdy*(dup + pav*du*gamma_minus_1)
                qzpo(i,j,kc,qpres ) = pnewrz
                qzpo(i,j,kc,qreint) = qzpo(i,j,kc,qpres) / gamma_minus_1
            end if

            ! ************************************************************************

            ! ************************************************************************
            rrlz = qzm(i,j,kc,qrho)
            rulz = rrlz*qzm(i,j,kc,qu)
            rvlz = rrlz*qzm(i,j,kc,qv)
            rwlz = rrlz*qzm(i,j,kc,qw)
            ekenlz = half*rrlz*(qzm(i,j,kc,qu)**2 + qzm(i,j,kc,qv)**2 &
                 + qzm(i,j,kc,qw)**2)
            relz = qzm(i,j,kc,qreint) + ekenlz

            ! Add transverse terms
            rrnewlz = rrlz - cdtdy*(fy(i,j+1,km,urho ) - fy(i,j,km,urho ))
            runewlz = rulz - cdtdy*(fy(i,j+1,km,umx  ) - fy(i,j,km,umx  ))
            rvnewlz = rvlz - cdtdy*(fy(i,j+1,km,umy  ) - fy(i,j,km,umy  ))
            rwnewlz = rwlz - cdtdy*(fy(i,j+1,km,umz  ) - fy(i,j,km,umz  ))
            renewlz = relz - cdtdy*(fy(i,j+1,km,ueden) - fy(i,j,km,ueden))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewlz .lt. zero) then
                 rrnewlz = rrlz
                 runewlz = rulz
                 rvnewlz = rvlz
                 rwnewlz = rwlz
                 renewlz = relz
            end if

            qzmo(i,j,kc,qrho) = rrnewlz
            qzmo(i,j,kc,qu) = runewlz/qzmo(i,j,kc,qrho)
            qzmo(i,j,kc,qv) = rvnewlz/qzmo(i,j,kc,qrho)
            qzmo(i,j,kc,qw) = rwnewlz/qzmo(i,j,kc,qrho)
            rhoekenlz = half*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,qrho)

            qzmo(i,j,kc,qreint)= renewlz - rhoekenlz
            qzmo(i,j,kc,qpres) = qzmo(i,j,kc,qreint) * gamma_minus_1

            if (qzmo(i,j,kc,qpres) .lt. small_pres) then
                pgp = pgdnvy(i,j+1,km)
                pgm = pgdnvy(i,j,km)
                ugp = ugdnvy(i,j+1,km)
                ugm = ugdnvy(i,j,km)
                dup = pgp*ugp - pgm*ugm
                pav = half*(pgp+pgm)
                du = ugp-ugm
                pnewlz = qzm(i,j,kc,qpres) - cdtdy*(dup + pav*du*gamma_minus_1)
                qzmo(i,j,kc,qpres ) = pnewlz
                qzmo(i,j,kc,qreint) = qzmo(i,j,kc,qpres) / gamma_minus_1
            end if
            ! ************************************************************************

         enddo
      enddo

      end subroutine transy2

      subroutine transz(qxm,qxmo,qxp,qxpo, &
                        qym,qymo,qyp,qypo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        fz,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                        ugdnvz,pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                        cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

      use amrex_mempool_module, only: amrex_allocate, amrex_deallocate
      use meth_params_module, only : qvar, nvar, qrho, qu, qv, qw, &
                                     qpres, qreint, &
                                     urho, umx, umy, umz, ueden, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
      integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      real(rt)  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) fz(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      real(rt) ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      real(rt) pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)

      ! Note that cdtdz = dtdz/3.d0/a_half
      real(rt) cdtdz

      integer n, nq
      integer i, j

      real(rt) rrnew, rr
      real(rt) compn, compu
      real(rt), pointer :: rrrx(:,:), rrry(:,:), rrlx(:,:), rrly(:,:)
      real(rt), pointer :: rurx(:,:), rury(:,:), rulx(:,:), ruly(:,:)
      real(rt), pointer :: rvrx(:,:), rvry(:,:), rvlx(:,:), rvly(:,:)
      real(rt), pointer :: rwrx(:,:), rwry(:,:), rwlx(:,:), rwly(:,:)
      real(rt), pointer :: ekenrx(:,:), ekenry(:,:), ekenlx(:,:), ekenly(:,:)
      real(rt), pointer :: rerx(:,:), rery(:,:), relx(:,:), rely(:,:)
      real(rt), pointer :: rrnewrx(:,:), rrnewry(:,:), rrnewlx(:,:), rrnewly(:,:)
      real(rt), pointer :: runewrx(:,:), runewry(:,:), runewlx(:,:), runewly(:,:)
      real(rt), pointer :: rvnewrx(:,:), rvnewry(:,:), rvnewlx(:,:), rvnewly(:,:)
      real(rt), pointer :: rwnewrx(:,:), rwnewry(:,:), rwnewlx(:,:), rwnewly(:,:)
      real(rt), pointer :: renewrx(:,:), renewry(:,:), renewlx(:,:), renewly(:,:)
      real(rt), pointer :: pnewrx(:,:),  pnewry(:,:),  pnewlx(:,:),  pnewly(:,:)
      real(rt), pointer :: rhoekenrx(:,:), rhoekenry(:,:), rhoekenlx(:,:), rhoekenly(:,:)
      real(rt), pointer :: pgp(:,:), pgm(:,:), ugp(:,:), ugm(:,:), dup(:,:), pav(:,:), du(:,:)

      integer          :: ipassive

      call amrex_allocate(rrrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rrry, ilo, ihi, jlo, jhi)
      call amrex_allocate(rrlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rrly, ilo, ihi, jlo, jhi)
      call amrex_allocate(rurx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rury, ilo, ihi, jlo, jhi)
      call amrex_allocate(rulx, ilo, ihi, jlo, jhi)
      call amrex_allocate(ruly, ilo, ihi, jlo, jhi)
      call amrex_allocate(rvrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rvry, ilo, ihi, jlo, jhi)
      call amrex_allocate(rvlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rvly, ilo, ihi, jlo, jhi)
      call amrex_allocate(rwrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rwry, ilo, ihi, jlo, jhi)
      call amrex_allocate(rwlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rwly, ilo, ihi, jlo, jhi)
      call amrex_allocate(ekenrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(ekenry, ilo, ihi, jlo, jhi)
      call amrex_allocate(ekenlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(ekenly, ilo, ihi, jlo, jhi)
      call amrex_allocate(rerx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rery, ilo, ihi, jlo, jhi)
      call amrex_allocate(relx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rely, ilo, ihi, jlo, jhi)
      call amrex_allocate(rrnewrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rrnewry, ilo, ihi, jlo, jhi)
      call amrex_allocate(rrnewlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rrnewly, ilo, ihi, jlo, jhi)
      call amrex_allocate(runewrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(runewry, ilo, ihi, jlo, jhi)
      call amrex_allocate(runewlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(runewly, ilo, ihi, jlo, jhi)
      call amrex_allocate(rvnewrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rvnewry, ilo, ihi, jlo, jhi)
      call amrex_allocate(rvnewlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rvnewly, ilo, ihi, jlo, jhi)
      call amrex_allocate(rwnewrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rwnewry, ilo, ihi, jlo, jhi)
      call amrex_allocate(rwnewlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rwnewly, ilo, ihi, jlo, jhi)
      call amrex_allocate(renewrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(renewry, ilo, ihi, jlo, jhi)
      call amrex_allocate(renewlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(renewly, ilo, ihi, jlo, jhi)
      call amrex_allocate(pnewrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(pnewry, ilo, ihi, jlo, jhi)
      call amrex_allocate(pnewlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(pnewly, ilo, ihi, jlo, jhi)
      call amrex_allocate(rhoekenrx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rhoekenry, ilo, ihi, jlo, jhi)
      call amrex_allocate(rhoekenlx, ilo, ihi, jlo, jhi)
      call amrex_allocate(rhoekenly, ilo, ihi, jlo, jhi)
      call amrex_allocate(pgp, ilo, ihi, jlo, jhi)
      call amrex_allocate(pgm, ilo, ihi, jlo, jhi)
      call amrex_allocate(ugp, ilo, ihi, jlo, jhi)
      call amrex_allocate(ugm, ilo, ihi, jlo, jhi)
      call amrex_allocate(dup, ilo, ihi, jlo, jhi)
      call amrex_allocate(pav, ilo, ihi, jlo, jhi)
      call amrex_allocate(du, ilo, ihi, jlo, jhi)

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            do i = ilo, ihi

                 compn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

                 if (i.ge.ilo+1) then
                    rr = qxp(i,j,km,qrho)
                    rrnew = rr - cdtdz*(fz(i,j,kc,urho) - fz(i,j,km,urho))
                    compu = rr*qxp(i,j,km,nq) - compn
                    qxpo(i,j,km,nq) = compu/rrnew
                 end if

                 if (j.ge.jlo+1) then
                    rr = qyp(i,j,km,qrho)
                    rrnew = rr - cdtdz*(fz(i,j,kc,urho) - fz(i,j,km,urho))
                    compu = rr*qyp(i,j,km,nq) - compn
                    qypo(i,j,km,nq) = compu/rrnew
                 end if

                 if (i.le.ihi-1) then
                    rr = qxm(i+1,j,km,qrho)
                    rrnew = rr - cdtdz*(fz(i,j,kc,urho) - fz(i,j,km,urho))
                    compu = rr*qxm(i+1,j,km,nq) - compn
                    qxmo(i+1,j,km,nq) = compu/rrnew
                 end if

                 if (j.le.jhi-1) then
                    rr = qym(i,j+1,km,qrho)
                    rrnew = rr - cdtdz*(fz(i,j,kc,urho) - fz(i,j,km,urho))
                    compu = rr*qym(i,j+1,km,nq) - compn
                    qymo(i,j+1,km,nq) = compu/rrnew
                 end if

            enddo
         enddo
      enddo

      pgp(ilo:ihi,jlo:jhi) = pgdnvz(ilo:ihi,jlo:jhi,kc)
      pgm(ilo:ihi,jlo:jhi) = pgdnvz(ilo:ihi,jlo:jhi,km)
      ugp(ilo:ihi,jlo:jhi) = ugdnvz(ilo:ihi,jlo:jhi,kc)
      ugm(ilo:ihi,jlo:jhi) = ugdnvz(ilo:ihi,jlo:jhi,km)


      ! Convert to conservation form
      rrrx(ilo+1:ihi,jlo:jhi) = qxp(ilo+1:ihi,jlo:jhi,km,qrho)
      rurx(ilo+1:ihi,jlo:jhi) = rrrx(ilo+1:ihi,jlo:jhi)*qxp(ilo+1:ihi,jlo:jhi,km,qu)
      rvrx(ilo+1:ihi,jlo:jhi) = rrrx(ilo+1:ihi,jlo:jhi)*qxp(ilo+1:ihi,jlo:jhi,km,qv)
      rwrx(ilo+1:ihi,jlo:jhi) = rrrx(ilo+1:ihi,jlo:jhi)*qxp(ilo+1:ihi,jlo:jhi,km,qw)
      ekenrx(ilo+1:ihi,jlo:jhi) = half*rrrx(ilo+1:ihi,jlo:jhi)*(qxp(ilo+1:ihi,jlo:jhi,km,qu)**2 + qxp(ilo+1:ihi,jlo:jhi,km,qv)**2 &
           + qxp(ilo+1:ihi,jlo:jhi,km,qw)**2)
      rerx(ilo+1:ihi,jlo:jhi) = qxp(ilo+1:ihi,jlo:jhi,km,qreint) + ekenrx(ilo+1:ihi,jlo:jhi)

      ! Add transverse terms
      rrnewrx(ilo+1:ihi,jlo:jhi) = rrrx(ilo+1:ihi,jlo:jhi) - cdtdz*(fz(ilo+1:ihi,jlo:jhi,kc,urho ) - fz(ilo+1:ihi,jlo:jhi,km,urho ))
      runewrx(ilo+1:ihi,jlo:jhi) = rurx(ilo+1:ihi,jlo:jhi) - cdtdz*(fz(ilo+1:ihi,jlo:jhi,kc,umx  ) - fz(ilo+1:ihi,jlo:jhi,km,umx  ))
      rvnewrx(ilo+1:ihi,jlo:jhi) = rvrx(ilo+1:ihi,jlo:jhi) - cdtdz*(fz(ilo+1:ihi,jlo:jhi,kc,umy  ) - fz(ilo+1:ihi,jlo:jhi,km,umy  ))
      rwnewrx(ilo+1:ihi,jlo:jhi) = rwrx(ilo+1:ihi,jlo:jhi) - cdtdz*(fz(ilo+1:ihi,jlo:jhi,kc,umz  ) - fz(ilo+1:ihi,jlo:jhi,km,umz  ))
      renewrx(ilo+1:ihi,jlo:jhi) = rerx(ilo+1:ihi,jlo:jhi) - cdtdz*(fz(ilo+1:ihi,jlo:jhi,kc,ueden) - fz(ilo+1:ihi,jlo:jhi,km,ueden))

      do j = jlo, jhi
          do i = ilo+1, ihi
             ! Reset to original value if adding transverse terms made density negative
             if (rrnewrx(i,j) .lt. zero) then
                rrnewrx(i,j) = rrrx(i,j)
                runewrx(i,j) = rurx(i,j)
                rvnewrx(i,j) = rvrx(i,j)
                rwnewrx(i,j) = rwrx(i,j)
                renewrx(i,j) = rerx(i,j)
             end if
          enddo
      enddo

      rrry(ilo:ihi,jlo+1:jhi) = qyp(ilo:ihi,jlo+1:jhi,km,qrho)
      rury(ilo:ihi,jlo+1:jhi) = rrry(ilo:ihi,jlo+1:jhi)*qyp(ilo:ihi,jlo+1:jhi,km,qu)
      rvry(ilo:ihi,jlo+1:jhi) = rrry(ilo:ihi,jlo+1:jhi)*qyp(ilo:ihi,jlo+1:jhi,km,qv)
      rwry(ilo:ihi,jlo+1:jhi) = rrry(ilo:ihi,jlo+1:jhi)*qyp(ilo:ihi,jlo+1:jhi,km,qw)
      ekenry(ilo:ihi,jlo+1:jhi) = half*rrry(ilo:ihi,jlo+1:jhi)*(qyp(ilo:ihi,jlo+1:jhi,km,qu)**2 + qyp(ilo:ihi,jlo+1:jhi,km,qv)**2 &
           + qyp(ilo:ihi,jlo+1:jhi,km,qw)**2)
      rery(ilo:ihi,jlo+1:jhi) = qyp(ilo:ihi,jlo+1:jhi,km,qreint) + ekenry(ilo:ihi,jlo+1:jhi)

      ! Add transverse terms
      rrnewry(ilo:ihi,jlo+1:jhi) = rrry(ilo:ihi,jlo+1:jhi) - cdtdz*(fz(ilo:ihi,jlo+1:jhi,kc,urho ) - fz(ilo:ihi,jlo+1:jhi,km,urho ))
      runewry(ilo:ihi,jlo+1:jhi) = rury(ilo:ihi,jlo+1:jhi) - cdtdz*(fz(ilo:ihi,jlo+1:jhi,kc,umx  ) - fz(ilo:ihi,jlo+1:jhi,km,umx  ))
      rvnewry(ilo:ihi,jlo+1:jhi) = rvry(ilo:ihi,jlo+1:jhi) - cdtdz*(fz(ilo:ihi,jlo+1:jhi,kc,umy  ) - fz(ilo:ihi,jlo+1:jhi,km,umy  ))
      rwnewry(ilo:ihi,jlo+1:jhi) = rwry(ilo:ihi,jlo+1:jhi) - cdtdz*(fz(ilo:ihi,jlo+1:jhi,kc,umz  ) - fz(ilo:ihi,jlo+1:jhi,km,umz  ))
      renewry(ilo:ihi,jlo+1:jhi) = rery(ilo:ihi,jlo+1:jhi) - cdtdz*(fz(ilo:ihi,jlo+1:jhi,kc,ueden) - fz(ilo:ihi,jlo+1:jhi,km,ueden))

      do j = jlo+1, jhi
          do i = ilo, ihi
              ! Reset to original value if adding transverse terms made density negative
              if (rrnewry(i,j) .lt. zero) then
                 rrnewry(i,j) = rrry(i,j)
                 runewry(i,j) = rury(i,j)
                 rvnewry(i,j) = rvry(i,j)
                 rwnewry(i,j) = rwry(i,j)
                 renewry(i,j) = rery(i,j)
              end if
          enddo
      enddo

      rrlx(ilo:ihi-1,jlo:jhi) = qxm(ilo+1:ihi,jlo:jhi,km,qrho)
      rulx(ilo:ihi-1,jlo:jhi) = rrlx(ilo:ihi-1,jlo:jhi)*qxm(ilo+1:ihi,jlo:jhi,km,qu)
      rvlx(ilo:ihi-1,jlo:jhi) = rrlx(ilo:ihi-1,jlo:jhi)*qxm(ilo+1:ihi,jlo:jhi,km,qv)
      rwlx(ilo:ihi-1,jlo:jhi) = rrlx(ilo:ihi-1,jlo:jhi)*qxm(ilo+1:ihi,jlo:jhi,km,qw)
      ekenlx(ilo:ihi-1,jlo:jhi) = half*rrlx(ilo:ihi-1,jlo:jhi)*(qxm(ilo+1:ihi,jlo:jhi,km,qu)**2 + qxm(ilo+1:ihi,jlo:jhi,km,qv)**2 &
           + qxm(ilo+1:ihi,jlo:jhi,km,qw)**2)
      relx(ilo:ihi-1,jlo:jhi) = qxm(ilo+1:ihi,jlo:jhi,km,qreint) + ekenlx(ilo:ihi-1,jlo:jhi)

      ! Add transverse terms
      rrnewlx(ilo:ihi-1,jlo:jhi) = rrlx(ilo:ihi-1,jlo:jhi) - cdtdz*(fz(ilo:ihi-1,jlo:jhi,kc,urho ) - fz(ilo:ihi-1,jlo:jhi,km,urho ))
      runewlx(ilo:ihi-1,jlo:jhi) = rulx(ilo:ihi-1,jlo:jhi) - cdtdz*(fz(ilo:ihi-1,jlo:jhi,kc,umx  ) - fz(ilo:ihi-1,jlo:jhi,km,umx  ))
      rvnewlx(ilo:ihi-1,jlo:jhi) = rvlx(ilo:ihi-1,jlo:jhi) - cdtdz*(fz(ilo:ihi-1,jlo:jhi,kc,umy  ) - fz(ilo:ihi-1,jlo:jhi,km,umy  ))
      rwnewlx(ilo:ihi-1,jlo:jhi) = rwlx(ilo:ihi-1,jlo:jhi) - cdtdz*(fz(ilo:ihi-1,jlo:jhi,kc,umz  ) - fz(ilo:ihi-1,jlo:jhi,km,umz  ))
      renewlx(ilo:ihi-1,jlo:jhi) = relx(ilo:ihi-1,jlo:jhi) - cdtdz*(fz(ilo:ihi-1,jlo:jhi,kc,ueden) - fz(ilo:ihi-1,jlo:jhi,km,ueden))

      do j = jlo, jhi
          do i = ilo, ihi-1
             ! Reset to original value if adding transverse terms made density negative
             if (rrnewlx(i,j) .lt. zero) then
                rrnewlx(i,j) = rrlx(i,j)
                runewlx(i,j) = rulx(i,j)
                rvnewlx(i,j) = rvlx(i,j)
                rwnewlx(i,j) = rwlx(i,j)
                renewlx(i,j) = relx(i,j)
             end if
          enddo
      enddo

      rrly(ilo:ihi,jlo:jhi-1) = qym(ilo:ihi,jlo+1:jhi,km,qrho)
      ruly(ilo:ihi,jlo:jhi-1) = rrly(ilo:ihi,jlo:jhi-1)*qym(ilo:ihi,jlo+1:jhi,km,qu)
      rvly(ilo:ihi,jlo:jhi-1) = rrly(ilo:ihi,jlo:jhi-1)*qym(ilo:ihi,jlo+1:jhi,km,qv)
      rwly(ilo:ihi,jlo:jhi-1) = rrly(ilo:ihi,jlo:jhi-1)*qym(ilo:ihi,jlo+1:jhi,km,qw)
      ekenly(ilo:ihi,jlo:jhi-1) = half*rrly(ilo:ihi,jlo:jhi-1)*(qym(ilo:ihi,jlo+1:jhi,km,qu)**2 + qym(ilo:ihi,jlo+1:jhi,km,qv)**2 &
           + qym(ilo:ihi,jlo+1:jhi,km,qw)**2)
      rely(ilo:ihi,jlo:jhi-1) = qym(ilo:ihi,jlo+1:jhi,km,qreint) + ekenly(ilo:ihi,jlo:jhi-1)

      ! Add transverse terms
      rrnewly(ilo:ihi,jlo:jhi-1) = rrly(ilo:ihi,jlo:jhi-1) - cdtdz*(fz(ilo:ihi,jlo:jhi-1,kc,urho ) - fz(ilo:ihi,jlo:jhi-1,km,urho ))
      runewly(ilo:ihi,jlo:jhi-1) = ruly(ilo:ihi,jlo:jhi-1) - cdtdz*(fz(ilo:ihi,jlo:jhi-1,kc,umx  ) - fz(ilo:ihi,jlo:jhi-1,km,umx  ))
      rvnewly(ilo:ihi,jlo:jhi-1) = rvly(ilo:ihi,jlo:jhi-1) - cdtdz*(fz(ilo:ihi,jlo:jhi-1,kc,umy  ) - fz(ilo:ihi,jlo:jhi-1,km,umy  ))
      rwnewly(ilo:ihi,jlo:jhi-1) = rwly(ilo:ihi,jlo:jhi-1) - cdtdz*(fz(ilo:ihi,jlo:jhi-1,kc,umz  ) - fz(ilo:ihi,jlo:jhi-1,km,umz  ))
      renewly(ilo:ihi,jlo:jhi-1) = rely(ilo:ihi,jlo:jhi-1) - cdtdz*(fz(ilo:ihi,jlo:jhi-1,kc,ueden) - fz(ilo:ihi,jlo:jhi-1,km,ueden))

      do j = jlo, jhi
          do i = ilo, ihi

              if (j.le.jhi-1) then
                ! Reset to original value if adding transverse terms made density negative
                if (rrnewly(i,j) .lt. zero) then
                   rrnewly(i,j) = rrly(i,j)
                   runewly(i,j) = ruly(i,j)
                   rvnewly(i,j) = rvly(i,j)
                   rwnewly(i,j) = rwly(i,j)
                   renewly(i,j) = rely(i,j)
                end if
              end if

          enddo
      enddo

      dup(ilo:ihi,jlo:jhi) = pgp(ilo:ihi,jlo:jhi)*ugp(ilo:ihi,jlo:jhi) - pgm(ilo:ihi,jlo:jhi)*ugm(ilo:ihi,jlo:jhi)
      pav(ilo:ihi,jlo:jhi) = half*(pgp(ilo:ihi,jlo:jhi)+pgm(ilo:ihi,jlo:jhi))
      du(ilo:ihi,jlo:jhi) = ugp(ilo:ihi,jlo:jhi)-ugm(ilo:ihi,jlo:jhi)


      ! Convert back to non-conservation form
      qxpo(ilo+1:ihi,jlo:jhi,km,qrho) = rrnewrx(ilo+1:ihi,jlo:jhi)
      qxpo(ilo+1:ihi,jlo:jhi,km,qu) = runewrx(ilo+1:ihi,jlo:jhi)/qxpo(ilo+1:ihi,jlo:jhi,km,qrho)
      qxpo(ilo+1:ihi,jlo:jhi,km,qv) = rvnewrx(ilo+1:ihi,jlo:jhi)/qxpo(ilo+1:ihi,jlo:jhi,km,qrho)
      qxpo(ilo+1:ihi,jlo:jhi,km,qw) = rwnewrx(ilo+1:ihi,jlo:jhi)/qxpo(ilo+1:ihi,jlo:jhi,km,qrho)
      rhoekenrx(ilo+1:ihi,jlo:jhi) = half*(runewrx(ilo+1:ihi,jlo:jhi)**2 &
                                         + rvnewrx(ilo+1:ihi,jlo:jhi)**2 &
                                         + rwnewrx(ilo+1:ihi,jlo:jhi)**2) / &
                                     qxpo(ilo+1:ihi,jlo:jhi,km,qrho)

      qxpo(ilo+1:ihi,jlo:jhi,km,qreint)= renewrx(ilo+1:ihi,jlo:jhi) - rhoekenrx(ilo+1:ihi,jlo:jhi)
      qxpo(ilo+1:ihi,jlo:jhi,km,qpres) = qxpo(ilo+1:ihi,jlo:jhi,km,qreint) * gamma_minus_1


      do j = jlo, jhi
          do i = ilo, ihi

              if (i.ge.ilo+1) then
                 if (qxpo(i,j,km,qpres) .lt. small_pres) then
                     pnewrx = qxp(i  ,j,km,qpres) - cdtdz*(dup(i,j) + pav(i,j)*du(i,j)*gamma_minus_1)
                     qxpo(i,j,km,qpres ) = pnewrx(i,j)
                     qxpo(i,j,km,qreint) = qxpo(i,j,km,qpres) / gamma_minus_1
                 end if
              end if

          enddo
      enddo

      qypo(ilo:ihi,jlo+1:jhi,km,qrho) = rrnewry(ilo:ihi,jlo+1:jhi)
      qypo(ilo:ihi,jlo+1:jhi,km,qu) = runewry(ilo:ihi,jlo+1:jhi)/qypo(ilo:ihi,jlo+1:jhi,km,qrho)
      qypo(ilo:ihi,jlo+1:jhi,km,qv) = rvnewry(ilo:ihi,jlo+1:jhi)/qypo(ilo:ihi,jlo+1:jhi,km,qrho)
      qypo(ilo:ihi,jlo+1:jhi,km,qw) = rwnewry(ilo:ihi,jlo+1:jhi)/qypo(ilo:ihi,jlo+1:jhi,km,qrho)
      rhoekenry(ilo:ihi,jlo+1:jhi) = half*(runewry(ilo:ihi,jlo+1:jhi)**2 &
                                         + rvnewry(ilo:ihi,jlo+1:jhi)**2 &
                                         + rwnewry(ilo:ihi,jlo+1:jhi)**2) / &
                                     qypo(ilo:ihi,jlo+1:jhi,km,qrho)

      qypo(ilo:ihi,jlo+1:jhi,km,qreint)= renewry(ilo:ihi,jlo+1:jhi) - rhoekenry(ilo:ihi,jlo+1:jhi)
      qypo(ilo:ihi,jlo+1:jhi,km,qpres) = qypo(ilo:ihi,jlo+1:jhi,km,qreint) * gamma_minus_1

      do j = jlo, jhi
          do i = ilo, ihi

              if (j.ge.jlo+1) then
                 if (qypo(i,j,km,qpres) .lt. small_pres) then
                     pnewry(i,j) = qyp(i,j  ,km,qpres) - cdtdz*(dup(i,j) + pav(i,j)*du(i,j)*gamma_minus_1)
                     qypo(i,j,km,qpres ) = pnewry(i,j)
                     qypo(i,j,km,qreint) = qypo(i,j,km,qpres) / gamma_minus_1
                 end if
              end if

          enddo
      enddo

      qxmo(ilo+1:ihi,jlo:jhi,km,qrho) = rrnewlx(ilo:ihi-1,jlo:jhi)
      qxmo(ilo+1:ihi,jlo:jhi,km,qu) = runewlx(ilo:ihi-1,jlo:jhi)/qxmo(ilo+1:ihi,jlo:jhi,km,qrho)
      qxmo(ilo+1:ihi,jlo:jhi,km,qv) = rvnewlx(ilo:ihi-1,jlo:jhi)/qxmo(ilo+1:ihi,jlo:jhi,km,qrho)
      qxmo(ilo+1:ihi,jlo:jhi,km,qw) = rwnewlx(ilo:ihi-1,jlo:jhi)/qxmo(ilo+1:ihi,jlo:jhi,km,qrho)
      rhoekenlx(ilo:ihi-1,jlo:jhi) = half*(runewlx(ilo:ihi-1,jlo:jhi)**2 &
                                         + rvnewlx(ilo:ihi-1,jlo:jhi)**2 &
                                         + rwnewlx(ilo:ihi-1,jlo:jhi)**2) / &
                                     qxmo(ilo+1:ihi,jlo:jhi,km,qrho)

      qxmo(ilo+1:ihi,jlo:jhi,km,qreint)= renewlx(ilo:ihi-1,jlo:jhi) - rhoekenlx(ilo:ihi-1,jlo:jhi)
      qxmo(ilo+1:ihi,jlo:jhi,km,qpres) = qxmo(ilo+1:ihi,jlo:jhi,km,qreint) * gamma_minus_1

      do j = jlo, jhi
          do i = ilo, ihi

              if (i.le.ihi-1) then
                 if (qxmo(i+1,j,km,qpres) .lt. small_pres) then
                     pnewlx(i,j) = qxm(i+1,j,km,qpres) - cdtdz*(dup(i,j) + pav(i,j)*du(i,j)*gamma_minus_1)
                     qxmo(i+1,j,km,qpres ) = pnewlx(i,j)
                     qxmo(i+1,j,km,qreint) = qxmo(i+1,j,km,qpres)  / gamma_minus_1
                 end if
              end if

          enddo
      enddo

      qymo(ilo:ihi,jlo+1:jhi,km,qrho) = rrnewly(ilo:ihi,jlo:jhi-1)
      qymo(ilo:ihi,jlo+1:jhi,km,qu) = runewly(ilo:ihi,jlo:jhi-1)/qymo(ilo:ihi,jlo+1:jhi,km,qrho)
      qymo(ilo:ihi,jlo+1:jhi,km,qv) = rvnewly(ilo:ihi,jlo:jhi-1)/qymo(ilo:ihi,jlo+1:jhi,km,qrho)
      qymo(ilo:ihi,jlo+1:jhi,km,qw) = rwnewly(ilo:ihi,jlo:jhi-1)/qymo(ilo:ihi,jlo+1:jhi,km,qrho)
      rhoekenly(ilo:ihi,jlo:jhi-1) = half*(runewly(ilo:ihi,jlo:jhi-1)**2 &
                                         + rvnewly(ilo:ihi,jlo:jhi-1)**2 &
                                         + rwnewly(ilo:ihi,jlo:jhi-1)**2) / &
                                     qymo(ilo:ihi,jlo+1:jhi,km,qrho)

      qymo(ilo:ihi,jlo+1:jhi,km,qreint)= renewly(ilo:ihi,jlo:jhi-1) - rhoekenly(ilo:ihi,jlo:jhi-1)
      qymo(ilo:ihi,jlo+1:jhi,km,qpres) = qymo(ilo:ihi,jlo+1:jhi,km,qreint) * gamma_minus_1

      do j = jlo, jhi
          do i = ilo, ihi

              if (j.le.jhi-1) then
                 if (qymo(i,j+1,km,qpres) .lt. small_pres) then
                     pnewly(i,j) = qym(i,j+1,km,qpres) - cdtdz*(dup(i,j) + pav(i,j)*du(i,j)*gamma_minus_1)
                     qymo(i,j+1,km,qpres ) = pnewly(i,j)
                     qymo(i,j+1,km,qreint) = qymo(i,j+1,km,qpres) / gamma_minus_1
                 endif
              end if

          enddo
      enddo

      call amrex_deallocate(rrrx)
      call amrex_deallocate(rrry)
      call amrex_deallocate(rrlx)
      call amrex_deallocate(rrly)
      call amrex_deallocate(rurx)
      call amrex_deallocate(rury)
      call amrex_deallocate(rulx)
      call amrex_deallocate(ruly)
      call amrex_deallocate(rvrx)
      call amrex_deallocate(rvry)
      call amrex_deallocate(rvlx)
      call amrex_deallocate(rvly)
      call amrex_deallocate(rwrx)
      call amrex_deallocate(rwry)
      call amrex_deallocate(rwlx)
      call amrex_deallocate(rwly)
      call amrex_deallocate(ekenrx)
      call amrex_deallocate(ekenry)
      call amrex_deallocate(ekenlx)
      call amrex_deallocate(ekenly)
      call amrex_deallocate(rerx)
      call amrex_deallocate(rery)
      call amrex_deallocate(relx)
      call amrex_deallocate(rely)
      call amrex_deallocate(rrnewrx)
      call amrex_deallocate(rrnewry)
      call amrex_deallocate(rrnewlx)
      call amrex_deallocate(rrnewly)
      call amrex_deallocate(runewrx)
      call amrex_deallocate(runewry)
      call amrex_deallocate(runewlx)
      call amrex_deallocate(runewly)
      call amrex_deallocate(rvnewrx)
      call amrex_deallocate(rvnewry)
      call amrex_deallocate(rvnewlx)
      call amrex_deallocate(rvnewly)
      call amrex_deallocate(rwnewrx)
      call amrex_deallocate(rwnewry)
      call amrex_deallocate(rwnewlx)
      call amrex_deallocate(rwnewly)
      call amrex_deallocate(renewrx)
      call amrex_deallocate(renewry)
      call amrex_deallocate(renewlx)
      call amrex_deallocate(renewly)
      call amrex_deallocate(pnewrx)
      call amrex_deallocate(pnewry)
      call amrex_deallocate(pnewlx)
      call amrex_deallocate(pnewly)
      call amrex_deallocate(rhoekenrx)
      call amrex_deallocate(rhoekenry)
      call amrex_deallocate(rhoekenlx)
      call amrex_deallocate(rhoekenly)
      call amrex_deallocate(pgp)
      call amrex_deallocate(pgm)
      call amrex_deallocate(ugp)
      call amrex_deallocate(ugm)
      call amrex_deallocate(dup)
      call amrex_deallocate(pav)
      call amrex_deallocate(du)

      end subroutine transz

      subroutine transxy(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fxy,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         fyx,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                         hdt,hdtdx,hdtdy,ilo,ihi,jlo,jhi,kc,km,k3d,a_old,a_new)

      use meth_params_module, only : qvar, nvar, qrho, qu, qv, qw, &
                                     qpres, qreint, &
                                     urho, umx, umy, umz, ueden, &
                                     small_pres, gamma_minus_1, &
                                     npassive, upass_map, qpass_map, & 
                                     version_2
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      real(rt)  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) fxy(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      real(rt) fyx(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      real(rt) ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      real(rt) pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      real(rt) ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      real(rt) pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      real(rt) srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      real(rt) hdt,a_old,a_new

      ! Note that hdtdx = dtdx/2.d0/a_half
      !       and hdtdy = dtdy/2.d0/a_half
      real(rt) hdtdx,hdtdy

      integer i, j
      integer n , nq

      real(rt) rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
      real(rt) rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
      real(rt) rrnewr, runewr, rvnewr, rwnewr, renewr, pnewr
      real(rt) rrnewl, runewl, rvnewl, rwnewl, renewl, pnewl
      real(rt) pgxp, pgxm, ugxp, ugxm
      real(rt) pgyp, pgym, ugyp, ugym
      real(rt) pgxpm, pgxmm, ugxpm, ugxmm
      real(rt) pgypm, pgymm, ugypm, ugymm
      real(rt) compr, compl, compnr, compnl
      real(rt) a_half

      real(rt) :: dux,duxm,duxp,duxpm,duy,duym,duyp,duypm
      real(rt) :: pxav,pxavm,pxnew,pxnewm,pyav,pyavm,pynew,pynewm
      integer          :: ipassive

      a_half = half * (a_old + a_new)

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            do i = ilo, ihi

               rrr = qp(i,j,kc,qrho)
               rrl = qm(i,j,kc,qrho)

               compr = rrr*qp(i,j,kc,nq)
               compl = rrl*qm(i,j,kc,nq)

               rrnewr = rrr - hdtdx*(fxy(i+1,j,kc,urho) - fxy(i,j,kc,urho)) &
                            - hdtdy*(fyx(i,j+1,kc,urho) - fyx(i,j,kc,urho))
               rrnewl = rrl - hdtdx*(fxy(i+1,j,km,urho) - fxy(i,j,km,urho)) &
                            - hdtdy*(fyx(i,j+1,km,urho) - fyx(i,j,km,urho))

               compnr = compr - hdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                              - hdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
               compnl = compl - hdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                              - hdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

               qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcq(i,j,k3d  ,nq) / a_half
               qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcq(i,j,k3d-1,nq) / a_half

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

             pgxp = pgdnvx(i+1,j,kc)
             pgxm = pgdnvx(i,j,kc)
             ugxp = ugdnvx(i+1,j,kc)
             ugxm = ugdnvx(i,j,kc)
 
             pgyp = pgdnvy(i,j+1,kc)
             pgym = pgdnvy(i,j,kc)
             ugyp = ugdnvy(i,j+1,kc)
             ugym = ugdnvy(i,j,kc)
 
             pgxpm = pgdnvx(i+1,j,km)
             pgxmm = pgdnvx(i,j,km)
             ugxpm = ugdnvx(i+1,j,km)
             ugxmm = ugdnvx(i,j,km)
 
             pgypm = pgdnvy(i,j+1,km)
             pgymm = pgdnvy(i,j,km)
             ugypm = ugdnvy(i,j+1,km)
             ugymm = ugdnvy(i,j,km)
 
             ! Convert to conservation form
             rrr = qp(i,j,kc,qrho)
             rur = rrr*qp(i,j,kc,qu)
             rvr = rrr*qp(i,j,kc,qv)
             rwr = rrr*qp(i,j,kc,qw)
             ekenr = half*rrr*(qp(i,j,kc,qu)**2 + qp(i,j,kc,qv)**2 + &
                  qp(i,j,kc,qw)**2)
             rer = qp(i,j,kc,qreint) + ekenr
 
             rrl = qm(i,j,kc,qrho)
             rul = rrl*qm(i,j,kc,qu)
             rvl = rrl*qm(i,j,kc,qv)
             rwl = rrl*qm(i,j,kc,qw)
             ekenl = half*rrl*(qm(i,j,kc,qu)**2 + qm(i,j,kc,qv)**2 + &
                  qm(i,j,kc,qw)**2)
             rel = qm(i,j,kc,qreint) + ekenl

             rrnewr = rrr - hdtdx*(fxy(i+1,j,kc,urho ) - fxy(i,j,kc,urho )) &
                          - hdtdy*(fyx(i,j+1,kc,urho ) - fyx(i,j,kc,urho ))
             runewr = rur - hdtdx*(fxy(i+1,j,kc,umx  ) - fxy(i,j,kc,umx  )) &
                          - hdtdy*(fyx(i,j+1,kc,umx  ) - fyx(i,j,kc,umx  ))
             rvnewr = rvr - hdtdx*(fxy(i+1,j,kc,umy  ) - fxy(i,j,kc,umy  )) &
                          - hdtdy*(fyx(i,j+1,kc,umy  ) - fyx(i,j,kc,umy  ))
             rwnewr = rwr - hdtdx*(fxy(i+1,j,kc,umz  ) - fxy(i,j,kc,umz  )) &
                          - hdtdy*(fyx(i,j+1,kc,umz  ) - fyx(i,j,kc,umz  ))
             renewr = rer - hdtdx*(fxy(i+1,j,kc,ueden) - fxy(i,j,kc,ueden)) &
                          - hdtdy*(fyx(i,j+1,kc,ueden) - fyx(i,j,kc,ueden))

             rrnewl = rrl - hdtdx*(fxy(i+1,j,km,urho ) - fxy(i,j,km,urho )) &
                          - hdtdy*(fyx(i,j+1,km,urho ) - fyx(i,j,km,urho ))
             runewl = rul - hdtdx*(fxy(i+1,j,km,umx  ) - fxy(i,j,km,umx  )) &
                          - hdtdy*(fyx(i,j+1,km,umx  ) - fyx(i,j,km,umx  ))
             rvnewl = rvl - hdtdx*(fxy(i+1,j,km,umy  ) - fxy(i,j,km,umy  )) &
                          - hdtdy*(fyx(i,j+1,km,umy  ) - fyx(i,j,km,umy  ))
             rwnewl = rwl - hdtdx*(fxy(i+1,j,km,umz  ) - fxy(i,j,km,umz  )) &
                          - hdtdy*(fyx(i,j+1,km,umz  ) - fyx(i,j,km,umz  ))
             renewl = rel - hdtdx*(fxy(i+1,j,km,ueden) - fxy(i,j,km,ueden)) &
                          - hdtdy*(fyx(i,j+1,km,ueden) - fyx(i,j,km,ueden))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewr .lt. zero) then
                 rrnewr = rrr
                 runewr = rur
                 rvnewr = rvr
                 rwnewr = rwr
                 renewr = rer
            end if
            if (rrnewl .lt. zero) then
                 rrnewl = rrl
                 runewl = rul
                 rvnewl = rvl
                 rwnewl = rwl
                 renewl = rel
            end if

            rhoekenr = half*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
            rhoekenl = half*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

            ! Convert back to non-conservation form
            ! ************************************************************************* 
            qpo(i,j,kc,qrho  ) = rrnewr

            qpo(i,j,kc,qu    ) = runewr/rrnewr
            qpo(i,j,kc,qv    ) = rvnewr/rrnewr
            qpo(i,j,kc,qw    ) = rwnewr/rrnewr

            qpo(i,j,kc,qreint) = renewr - rhoekenr + hdt* srcq(i,j,k3d,qreint) / a_old
            qpo(i,j,kc,qpres) = qpo(i,j,kc,qreint) * gamma_minus_1

            if (qpo(i,j,kc,qpres) .lt. small_pres) then
                duxp = pgxp*ugxp - pgxm*ugxm
                pxav = half*(pgxp+pgxm)
                dux = ugxp-ugxm
                pxnew = hdtdx*(duxp + pxav*dux*gamma_minus_1)

                duyp = pgyp*ugyp - pgym*ugym
                pyav = half*(pgyp+pgym)
                duy = ugyp-ugym
                pynew = hdtdy*(duyp + pyav*duy*gamma_minus_1)

                pnewr = qp(i,j,kc,qpres) - pxnew - pynew
                qpo(i,j,kc,qpres ) = pnewr + hdt* srcq(i,j,k3d,qpres ) / a_old
                qpo(i,j,kc,qreint) = qpo(i,j,kc,qpres) / gamma_minus_1
            end if

            ! ************************************************************************* 
            qmo(i,j,kc,qrho  ) = rrnewl

            qmo(i,j,kc,qu    ) = runewl/rrnewl
            qmo(i,j,kc,qv    ) = rvnewl/rrnewl
            qmo(i,j,kc,qw    ) = rwnewl/rrnewl

            qmo(i,j,kc,qreint) = renewl - rhoekenl + hdt* srcq(i,j,k3d-1,qreint) / a_old
            qmo(i,j,kc,qpres) = qmo(i,j,kc,qreint) * gamma_minus_1

            if (qmo(i,j,kc,qpres) .lt. small_pres) then

                duxpm = pgxpm*ugxpm - pgxmm*ugxmm
                pxavm = half*(pgxpm+pgxmm)
                duxm = ugxpm-ugxmm
                pxnewm = hdtdx*(duxpm + pxavm*duxm*gamma_minus_1)

                duypm = pgypm*ugypm - pgymm*ugymm
                pyavm = half*(pgypm+pgymm)
                duym = ugypm-ugymm
                pynewm = hdtdy*(duypm + pyavm*duym*gamma_minus_1)

                pnewl = qm(i,j,kc,qpres) - pxnewm - pynewm
                qmo(i,j,kc,qpres ) = pnewl + hdt* srcq(i,j,k3d-1,qpres ) / a_old
                qmo(i,j,kc,qreint) = qmo(i,j,kc,qpres) / gamma_minus_1
            end if
            ! ************************************************************************* 
         enddo
      enddo

      ! Version_2 = 0: we add the full source term here
      ! Version_2 = 1: we already added the (piecewise constant) source term in the original tracing
      ! Version_2 = 2: we already added the (parabolic         ) source term in the original tracing
      ! Version_2 = 3: we project then add the source term in trace*_src_3d

      if (version_2 .eq. 0) then
         do j = jlo, jhi
            do i = ilo, ihi
               qpo(i,j,kc,qrho) = qpo(i,j,kc,qrho) + hdt*srcq(i,j,k3d,qrho) / a_old
               qpo(i,j,kc,qu  ) = qpo(i,j,kc,qu  ) + hdt*srcq(i,j,k3d,qu  ) / a_old
               qpo(i,j,kc,qv  ) = qpo(i,j,kc,qv  ) + hdt*srcq(i,j,k3d,qv  ) / a_old
               qpo(i,j,kc,qw  ) = qpo(i,j,kc,qw  ) + hdt*srcq(i,j,k3d,qw  ) / a_old

               qmo(i,j,kc,qrho) = qmo(i,j,kc,qrho) + hdt*srcq(i,j,k3d-1,qrho) / a_old
               qmo(i,j,kc,qu  ) = qmo(i,j,kc,qu  ) + hdt*srcq(i,j,k3d-1,qu  ) / a_old
               qmo(i,j,kc,qv  ) = qmo(i,j,kc,qv  ) + hdt*srcq(i,j,k3d-1,qv  ) / a_old 
               qmo(i,j,kc,qw  ) = qmo(i,j,kc,qw  ) + hdt*srcq(i,j,k3d-1,qw  ) / a_old
            enddo
         enddo
      endif

      end subroutine transxy

      subroutine transxz(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fxz,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         fzx,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         ugdnvz,pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                         hdt,hdtdx,hdtdz,ilo,ihi,jlo,jhi,km,kc,k3d,a_old,a_new)

      use meth_params_module, only : qvar, nvar, qrho, qu, qv, qw, &
                                     qpres, qreint, &
                                     urho, umx, umy, umz, ueden, &
                                     small_pres, gamma_minus_1, &
                                     npassive, upass_map, qpass_map, & 
                                     version_2
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      real(rt)  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) fxz(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      real(rt) fzx(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      real(rt) ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      real(rt) pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      real(rt) ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      real(rt) pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      real(rt) srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      real(rt) hdt,a_old,a_new

      ! Note that hdtdx = dtdx/2.d0/a_half
      !       and hdtdz = dtdz/2.d0/a_half
      real(rt) hdtdx,hdtdz

      integer i, j
      integer n, nq

      real(rt) rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
      real(rt) rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
      real(rt) rrnewr, runewr, rvnewr, rwnewr, renewr, pnewr
      real(rt) rrnewl, runewl, rvnewl, rwnewl, renewl, pnewl
      real(rt) pgxp, pgxm, ugxp, ugxm
      real(rt) pgzp, pgzm, ugzp, ugzm
      real(rt) compr, compl, compnr, compnl
      real(rt) a_half

      real(rt) :: drr, dcompn
      real(rt) :: dux, duxp, duz, duzp, pxav, pxnew, pzav, pznew
      integer          :: ipassive

      a_half = half * (a_old + a_new)
 
      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            do i = ilo, ihi

               drr    = - hdtdx*(fxz(i+1,j,km,urho) - fxz(i,j,km,urho)) &
                        - hdtdz*(fzx(i  ,j,kc,urho) - fzx(i,j,km,urho))
               dcompn = - hdtdx*(fxz(i+1,j,km,n   ) - fxz(i,j,km,n)) &
                        - hdtdz*(fzx(i  ,j,kc,n   ) - fzx(i,j,km,n))

               if (j.ge.jlo+1) then
                  rrr = qp(i,j,km,qrho)
                  compr = rrr*qp(i,j,km,nq)

                  rrnewr = rrr   + drr
                  compnr = compr + dcompn

                  qpo(i,j  ,km,nq) = compnr/rrnewr + hdt*srcq(i,j,k3d,nq) / a_half
               end if

               if (j.le.jhi-1) then
                  rrl = qm(i,j+1,km,qrho)
                  compl = rrl*qm(i,j+1,km,nq)

                  rrnewl = rrl   + drr
                  compnl = compl + dcompn

                  qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcq(i,j,k3d,nq) / a_half
               end if

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

            pgxp = pgdnvx(i+1,j,km)
            pgxm = pgdnvx(i,j,km)
            ugxp = ugdnvx(i+1,j,km)
            ugxm = ugdnvx(i,j,km)

            pgzp = pgdnvz(i,j,kc)
            pgzm = pgdnvz(i,j,km)
            ugzp = ugdnvz(i,j,kc)
            ugzm = ugdnvz(i,j,km)

            if (j.ge.jlo+1) then
               ! Convert to conservation form
               rrr = qp(i,j,km,qrho)
               rur = rrr*qp(i,j,km,qu)
               rvr = rrr*qp(i,j,km,qv)
               rwr = rrr*qp(i,j,km,qw)
               ekenr = half*rrr*(qp(i,j,km,qu)**2 + qp(i,j,km,qv)**2 + qp(i,j,km,qw)**2)
               rer = qp(i,j,km,qreint) + ekenr

               ! Add transverse terms
               rrnewr = rrr - hdtdx*(fxz(i+1,j,km,urho ) - fxz(i,j,km,urho )) &
                            - hdtdz*(fzx(i  ,j,kc,urho ) - fzx(i,j,km,urho ))
               runewr = rur - hdtdx*(fxz(i+1,j,km,umx  ) - fxz(i,j,km,umx  )) &
                            - hdtdz*(fzx(i  ,j,kc,umx  ) - fzx(i,j,km,umx  ))
               rvnewr = rvr - hdtdx*(fxz(i+1,j,km,umy  ) - fxz(i,j,km,umy  )) &
                            - hdtdz*(fzx(i  ,j,kc,umy  ) - fzx(i,j,km,umy  ))
               rwnewr = rwr - hdtdx*(fxz(i+1,j,km,umz  ) - fxz(i,j,km,umz  )) &
                            - hdtdz*(fzx(i  ,j,kc,umz  ) - fzx(i,j,km,umz  ))
               renewr = rer - hdtdx*(fxz(i+1,j,km,ueden) - fxz(i,j,km,ueden)) &
                            - hdtdz*(fzx(i  ,j,kc,ueden) - fzx(i,j,km,ueden))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewr .lt. zero) then
                    rrnewr = rrr
                    runewr = rur
                    rvnewr = rvr
                    rwnewr = rwr
                    renewr = rer
               end if

               rhoekenr = half*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr

            end if

            if (j.le.jhi-1) then
               rrl = qm(i,j+1,km,qrho)
               rul = rrl*qm(i,j+1,km,qu)
               rvl = rrl*qm(i,j+1,km,qv)
               rwl = rrl*qm(i,j+1,km,qw)
               ekenl = half*rrl*(qm(i,j+1,km,qu)**2 + qm(i,j+1,km,qv)**2 + qm(i,j+1,km,qw)**2)
               rel = qm(i,j+1,km,qreint) + ekenl

              ! Add transverse terms
              rrnewl = rrl - hdtdx*(fxz(i+1,j,km,urho ) - fxz(i,j,km,urho )) &
                           - hdtdz*(fzx(i  ,j,kc,urho ) - fzx(i,j,km,urho ))
              runewl = rul - hdtdx*(fxz(i+1,j,km,umx  ) - fxz(i,j,km,umx  )) &
                           - hdtdz*(fzx(i  ,j,kc,umx  ) - fzx(i,j,km,umx  ))
              rvnewl = rvl - hdtdx*(fxz(i+1,j,km,umy  ) - fxz(i,j,km,umy  )) &
                           - hdtdz*(fzx(i  ,j,kc,umy  ) - fzx(i,j,km,umy  ))
              rwnewl = rwl - hdtdx*(fxz(i+1,j,km,umz  ) - fxz(i,j,km,umz  )) &
                           - hdtdz*(fzx(i  ,j,kc,umz  ) - fzx(i,j,km,umz  ))
              renewl = rel - hdtdx*(fxz(i+1,j,km,ueden) - fxz(i,j,km,ueden)) &
                           - hdtdz*(fzx(i  ,j,kc,ueden) - fzx(i,j,km,ueden))
               ! Reset to original value if adding transverse terms made density negative
               if (rrnewl .lt. zero) then
                    rrnewl = rrl
                    runewl = rul
                    rvnewl = rvl
                    rwnewl = rwl
                    renewl = rel
               end if

               rhoekenl = half*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

            end if

            ! These are used only if the (rho e) version goes negative 
            duxp = pgxp*ugxp - pgxm*ugxm
            pxav = half*(pgxp+pgxm)
            dux = ugxp-ugxm
            pxnew = hdtdx*(duxp + pxav*dux*gamma_minus_1)

            duzp = pgzp*ugzp - pgzm*ugzm
            pzav = half*(pgzp+pgzm)
            duz = ugzp-ugzm
            pznew = hdtdz*(duzp + pzav*duz*gamma_minus_1)
            ! End comment

            ! Convert back to non-conservation form
            ! ************************************************************************* 
            if (j.ge.jlo+1) then
               qpo(i,j,km,qrho  ) = rrnewr
               qpo(i,j,km,qu    ) = runewr/rrnewr
               qpo(i,j,km,qv    ) = rvnewr/rrnewr
               qpo(i,j,km,qw    ) = rwnewr/rrnewr

               qpo(i,j,km,qreint) = renewr - rhoekenr + hdt* srcq(i,j,k3d,qreint) / a_old
               qpo(i,j,km,qpres) = qpo(i,j,km,qreint) * gamma_minus_1

               if (qpo(i,j,km,qpres) .lt. small_pres) then
                   pnewr = qp(i,j  ,km,qpres) - pxnew - pznew
                   qpo(i,j,km,qpres ) = pnewr + hdt* srcq(i,j,k3d,qpres ) / a_old
                   qpo(i,j,km,qreint) = qpo(i,j,km,qpres) / gamma_minus_1
               end if
            end if

            ! ************************************************************************* 
            if (j.le.jhi-1) then
               qmo(i,j+1,km,qrho  ) = rrnewl
               qmo(i,j+1,km,qu    ) = runewl/rrnewl
               qmo(i,j+1,km,qv    ) = rvnewl/rrnewl
               qmo(i,j+1,km,qw    ) = rwnewl/rrnewl

               qmo(i,j+1,km,qreint) = renewl - rhoekenl + hdt* srcq(i,j,k3d,qreint) / a_old
               qmo(i,j+1,km,qpres) = qmo(i,j+1,km,qreint) * gamma_minus_1

               if (qmo(i,j+1,km,qpres) .lt. small_pres) then
                   pnewl = qm(i,j+1,km,qpres) - pxnew - pznew
                   qmo(i,j+1,km,qpres ) = pnewl + hdt* srcq(i,j,k3d,qpres ) / a_old
                   qmo(i,j+1,km,qreint) = qmo(i,j+1,km,qpres) / gamma_minus_1
               end if

            end if
            ! ************************************************************************* 

         enddo
      enddo

      ! Version_2 = 0: we add the full source term here
      ! Version_2 = 1: we already added the (piecewise constant) source term in the original tracing
      ! Version_2 = 2: we already added the (parabolic         ) source term in the original tracing
      ! Version_2 = 3: we project then add the source term in trace*_src_3d

      if (version_2 .eq. 0) then
         do j = jlo, jhi
            do i = ilo, ihi
               if (j.ge.jlo+1) then
                  qpo(i,j,km,qrho) = qpo(i,j,km,qrho) + hdt*srcq(i,j,k3d,qrho  ) / a_old
                  qpo(i,j,km,qu  ) = qpo(i,j,km,qu  ) + hdt*srcq(i,j,k3d,qu) / a_old
                  qpo(i,j,km,qv  ) = qpo(i,j,km,qv  ) + hdt*srcq(i,j,k3d,qv) / a_old
                  qpo(i,j,km,qw  ) = qpo(i,j,km,qw  ) + hdt*srcq(i,j,k3d,qw) / a_old
               end if

               if (j.le.jhi-1) then
                  qmo(i,j+1,km,qrho) = qmo(i,j+1,km,qrho) + hdt*srcq(i,j,k3d,qrho  ) / a_old
                  qmo(i,j+1,km,qu  ) = qmo(i,j+1,km,qu  ) + hdt*srcq(i,j,k3d,qu) / a_old
                  qmo(i,j+1,km,qv  ) = qmo(i,j+1,km,qv  ) + hdt*srcq(i,j,k3d,qv) / a_old
                  qmo(i,j+1,km,qw  ) = qmo(i,j+1,km,qw  ) + hdt*srcq(i,j,k3d,qw) / a_old
               end if
            enddo
         enddo
      endif

      end subroutine transxz

      subroutine transyz(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fyz,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         fzy,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         ugdnvz,pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                         hdt,hdtdy,hdtdz,ilo,ihi,jlo,jhi,km,kc,k3d,a_old,a_new)

      use meth_params_module, only : qvar, nvar, qrho, qu, qv, qw, &
                                     qpres, qreint, &
                                     urho, umx, umy, umz, ueden, &
                                     small_pres, gamma_minus_1, &
                                     npassive, upass_map, qpass_map, & 
                                     version_2
      implicit none

      integer :: qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer :: fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer :: fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
      integer :: pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer :: pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
      integer :: src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer :: ilo,ihi,jlo,jhi,km,kc,k3d

      real(rt) :: qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) :: qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) :: qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) :: qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt) :: fyz(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      real(rt) :: fzy(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      real(rt) :: ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      real(rt) :: pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      real(rt) :: ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      real(rt) :: pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      real(rt) :: srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      real(rt) :: hdt,a_old,a_new

      ! Note that hdtdy = dtdy/2.d0/a_half
      !       and hdtdz = dtdz/2.d0/a_half
      real(rt) :: hdtdy,hdtdz

      integer :: i, j
      integer :: n, nq

      real(rt) :: rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
      real(rt) :: rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
      real(rt) :: rrnewr, runewr, rvnewr, rwnewr, renewr, pnewr
      real(rt) :: rrnewl, runewl, rvnewl, rwnewl, renewl, pnewl
      real(rt) :: pgyp, pgym, ugyp, ugym
      real(rt) :: pgzp, pgzm, ugzp, ugzm
      real(rt) :: compr, compl, compnr, compnl
      real(rt) :: a_half

      real(rt) :: drr, dcompn
      real(rt) :: duy,duyp,duz,duzp,pyav,pynew,pzav,pznew
      integer          :: ipassive

      a_half = half * (a_old + a_new)

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            do i = ilo, ihi

               drr    = - hdtdy*(fyz(i,j+1,km,urho) - fyz(i,j,km,urho)) &
                        - hdtdz*(fzy(i,j  ,kc,urho) - fzy(i,j,km,urho))
               dcompn = - hdtdy*(fyz(i,j+1,km,n   ) - fyz(i,j,km,n)) &
                        - hdtdz*(fzy(i,j  ,kc,n   ) - fzy(i,j,km,n))

               if (i.ge.ilo+1) then
                  rrr = qp(i,j,km,qrho)
                  compr  = rrr*qp(i,j,km,nq)
   
                  rrnewr = rrr   + drr
                  compnr = compr + dcompn

                  qpo(i  ,j,km,nq) = compnr/rrnewr + hdt*srcq(i,j,k3d,nq) / a_half

               end if

               if (i.le.ihi-1) then
                  rrl = qm(i+1,j,km,qrho)
                  compl  = rrl*qm(i+1,j,km,nq)

                  rrnewl = rrl + drr
                  compnl = compl + dcompn

                  qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcq(i,j,k3d,nq) / a_half
               end if

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

            pgyp = pgdnvy(i,j+1,km)
            pgym = pgdnvy(i,j,km)
            ugyp = ugdnvy(i,j+1,km)
            ugym = ugdnvy(i,j,km)

            pgzp = pgdnvz(i,j,kc)
            pgzm = pgdnvz(i,j,km)
            ugzp = ugdnvz(i,j,kc)
            ugzm = ugdnvz(i,j,km)

            if (i.ge.ilo+1) then
               ! Convert to conservation form
               rrr =     qp(i,j,km,qrho)
               rur = rrr*qp(i,j,km,qu)
               rvr = rrr*qp(i,j,km,qv)
               rwr = rrr*qp(i,j,km,qw)
               ekenr = half*rrr*(qp(i,j,km,qu)**2 + qp(i,j,km,qv)**2 + &
                                  qp(i,j,km,qw)**2)
               rer = qp(i,j,km,qreint) + ekenr

               ! Add transverse terms
               rrnewr = rrr - hdtdy*(fyz(i,j+1,km,urho ) - fyz(i,j,km,urho )) &
                            - hdtdz*(fzy(i,j  ,kc,urho ) - fzy(i,j,km,urho ))
               runewr = rur - hdtdy*(fyz(i,j+1,km,umx  ) - fyz(i,j,km,umx  )) &
                            - hdtdz*(fzy(i,j  ,kc,umx  ) - fzy(i,j,km,umx  ))
               rvnewr = rvr - hdtdy*(fyz(i,j+1,km,umy  ) - fyz(i,j,km,umy  )) &
                            - hdtdz*(fzy(i,j  ,kc,umy  ) - fzy(i,j,km,umy  ))
               rwnewr = rwr - hdtdy*(fyz(i,j+1,km,umz  ) - fyz(i,j,km,umz  )) &
                            - hdtdz*(fzy(i,j  ,kc,umz  ) - fzy(i,j,km,umz  ))
               renewr = rer - hdtdy*(fyz(i,j+1,km,ueden) - fyz(i,j,km,ueden)) &
                            - hdtdz*(fzy(i,j  ,kc,ueden) - fzy(i,j,km,ueden))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewr .lt. zero) then
                    rrnewr = rrr
                    runewr = rur
                    rvnewr = rvr
                    rwnewr = rwr
                    renewr = rer
               end if

               rhoekenr = half*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr

            end if

            if (i.le.ihi-1) then
               rrl =     qm(i+1,j,km,qrho)
               rul = rrl*qm(i+1,j,km,qu)
               rvl = rrl*qm(i+1,j,km,qv)
               rwl = rrl*qm(i+1,j,km,qw)
               ekenl = half*rrl*(qm(i+1,j,km,qu)**2 + qm(i+1,j,km,qv)**2 + &
                                  qm(i+1,j,km,qw)**2)
               rel = qm(i+1,j,km,qreint) + ekenl

            ! Add transverse terms
               rrnewl = rrl - hdtdy*(fyz(i,j+1,km,urho ) - fyz(i,j,km,urho )) &
                            - hdtdz*(fzy(i,j  ,kc,urho ) - fzy(i,j,km,urho ))
               runewl = rul - hdtdy*(fyz(i,j+1,km,umx  ) - fyz(i,j,km,umx  )) &
                            - hdtdz*(fzy(i,j  ,kc,umx  ) - fzy(i,j,km,umx  ))
               rvnewl = rvl - hdtdy*(fyz(i,j+1,km,umy  ) - fyz(i,j,km,umy  )) &
                            - hdtdz*(fzy(i,j  ,kc,umy  ) - fzy(i,j,km,umy  ))
               rwnewl = rwl - hdtdy*(fyz(i,j+1,km,umz  ) - fyz(i,j,km,umz  )) &
                            - hdtdz*(fzy(i,j  ,kc,umz  ) - fzy(i,j,km,umz  ))
               renewl = rel - hdtdy*(fyz(i,j+1,km,ueden) - fyz(i,j,km,ueden)) &
                         - hdtdz*(fzy(i,j  ,kc,ueden) - fzy(i,j,km,ueden))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewl .lt. zero) then
                    rrnewl = rrl
                    runewl = rul
                    rvnewl = rvl
                    rwnewl = rwl
                    renewl = rel
               end if

               rhoekenl = half*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

            end if

            ! These are used only if the (rho e) version goes negative 
            duyp = pgyp*ugyp - pgym*ugym
            pyav = half*(pgyp+pgym)
            duy = ugyp-ugym
            pynew = hdtdy*(duyp + pyav*duy*gamma_minus_1)

            duzp = pgzp*ugzp - pgzm*ugzm
            pzav = half*(pgzp+pgzm)
            duz = ugzp-ugzm
            pznew = hdtdz*(duzp + pzav*duz*gamma_minus_1)
            ! End section

            ! Convert back to non-conservation form
            ! ************************************************************************* 
            if (i.ge.ilo+1) then
               qpo(i,j,km,qrho  ) = rrnewr
               qpo(i,j,km,qu    ) = runewr/rrnewr
               qpo(i,j,km,qv    ) = rvnewr/rrnewr
               qpo(i,j,km,qw    ) = rwnewr/rrnewr

               qpo(i,j,km,qreint) = renewr - rhoekenr + hdt* srcq(i,j,k3d,qreint) / a_old
               qpo(i,j,km,qpres) = qpo(i,j,km,qreint) * gamma_minus_1

               if (qpo(i,j,km,qpres) .lt. small_pres) then
                   pnewr = qp(i,j,km,qpres) - pynew - pznew
                   qpo(i,j,km,qpres ) = pnewr + hdt* srcq(i,j,k3d,qpres ) / a_old
                   qpo(i,j,km,qreint) = qpo(i,j,km,qpres) / gamma_minus_1
               end if
            end if

            ! ************************************************************************* 
            if (i.le.ihi-1) then
               qmo(i+1,j,km,qrho  ) = rrnewl
               qmo(i+1,j,km,qu    ) = runewl/rrnewl
               qmo(i+1,j,km,qv    ) = rvnewl/rrnewl
               qmo(i+1,j,km,qw    ) = rwnewl/rrnewl

               qmo(i+1,j,km,qreint) = renewl - rhoekenl + hdt* srcq(i,j,k3d,qreint) / a_old
               qmo(i+1,j,km,qpres) = qmo(i+1,j,km,qreint) * gamma_minus_1

               if (qmo(i+1,j,km,qpres) .lt. small_pres) then
                   pnewl = qm(i+1,j,km,qpres) - pynew - pznew
                   qmo(i+1,j,km,qpres ) = pnewl + hdt* srcq(i,j,k3d,qpres ) / a_old
                   qmo(i+1,j,km,qreint) = qmo(i+1,j,km,qpres) / gamma_minus_1
               end if
            end if
            ! ************************************************************************* 

         enddo
      enddo

      ! Version_2 = 0: we add the full source term here
      ! Version_2 = 1: we already added the (piecewise constant) source term in the original tracing
      ! Version_2 = 2: we already added the (parabolic         ) source term in the original tracing
      ! Version_2 = 3: we project then add the source term in trace*_src_3d

      if (version_2 .eq. 0) then
         do j = jlo, jhi
            do i = ilo, ihi
               if (i.ge.ilo+1) then
                  qpo(i,j,km,qrho) = qpo(i,j,km,qrho) + hdt*srcq(i,j,k3d,qrho  ) / a_old
                  qpo(i,j,km,qu  ) = qpo(i,j,km,qu  ) + hdt*srcq(i,j,k3d,qu) / a_old
                  qpo(i,j,km,qv  ) = qpo(i,j,km,qv  ) + hdt*srcq(i,j,k3d,qv) / a_old
                  qpo(i,j,km,qw  ) = qpo(i,j,km,qw  ) + hdt*srcq(i,j,k3d,qw) / a_old
               end if
   
               if (i.le.ihi-1) then
                  qmo(i+1,j,km,qrho) = qmo(i+1,j,km,qrho) + hdt*srcq(i,j,k3d,qrho  ) / a_old
                  qmo(i+1,j,km,qu  ) = qmo(i+1,j,km,qu  ) + hdt*srcq(i,j,k3d,qu) / a_old
                  qmo(i+1,j,km,qv  ) = qmo(i+1,j,km,qv  ) + hdt*srcq(i,j,k3d,qv) / a_old
                  qmo(i+1,j,km,qw  ) = qmo(i+1,j,km,qw  ) + hdt*srcq(i,j,k3d,qw) / a_old
               end if
            enddo
         enddo
      endif

      end subroutine transyz

end module transverse_module
````

