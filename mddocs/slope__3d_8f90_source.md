
# File slope\_3d.f90

[**File List**](files.md) **>** [**HydroFortran**](dir_1fab266cd447ad3f3624320661f845f1.md) **>** [**slope\_3d.f90**](slope__3d_8f90.md)

[Go to the documentation of this file.](slope__3d_8f90.md) 


````cpp
module slope_module

  implicit none

  private

  public uslope

contains

      subroutine uslope(q,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        dqx,dqy,dqz,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,nv)

      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module
      use amrex_constants_module

      implicit none

      integer ilo,ihi
      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2,kc,k3d,nv

      real(rt) q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,nv)
      real(rt) flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      real(rt) dqx(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)
      real(rt) dqy(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)
      real(rt) dqz(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)

      integer i, j, k, n

      real(rt) dlft, drgt, slop, dq1
      real(rt) dm, dp, dc, ds, sl, dl, dfm, dfp

      real(rt), allocatable::dsgn(:,:),dlim(:,:),df(:,:),dcen(:,:)

      ilo = min(ilo1,ilo2)
      ihi = max(ihi1,ihi2)

      allocate (dsgn(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2,ilo-2:ihi+2))

      if (iorder.eq.1) then

         do n = 1, nv
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqx(i,j,kc,n) = zero
                  dqy(i,j,kc,n) = zero
                  dqz(i,j,kc,n) = zero
               enddo
            enddo
         enddo

      else

         do n = 1, nv

            ! Compute slopes in first coordinate direction
            do j = ilo2-1, ihi2+1

               ! First compute Fromm slopes
               do i = ilo1-2, ihi1+2
                  dlft = two*(q(i ,j,k3d,n) - q(i-1,j,k3d,n))
                  drgt = two*(q(i+1,j,k3d,n) - q(i ,j,k3d,n))
                  dcen(i,j) = fourth * (dlft+drgt)
                  dsgn(i,j) = sign(one, dcen(i,j))
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. zero) then
                     dlim(i,j) = slop
                  else
                     dlim(i,j) = zero
                  endif
                  df(i,j) = dsgn(i,j)*min( dlim(i,j), abs(dcen(i,j)) )
               enddo

               ! Now compute limited fourth order slopes
               do i = ilo1-1, ihi1+1
                  dq1       = four3rd*dcen(i,j) - sixth*(df(i+1,j) + df(i-1,j))
                  dqx(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
               enddo

            enddo

            ! Compute slopes in second coordinate direction
            do i = ilo1-1, ihi1+1
               ! First compute Fromm slopes for this column
               do j = ilo2-2, ihi2+2
                  dlft = two*(q(i,j ,k3d,n) - q(i,j-1,k3d,n))
                  drgt = two*(q(i,j+1,k3d,n) - q(i,j ,k3d,n))
                  dcen(i,j) = fourth * (dlft+drgt)
                  dsgn(i,j) = sign( one, dcen(i,j) )
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. zero) then
                     dlim(i,j) = slop
                  else
                     dlim(i,j) = zero
                  endif
                  df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
               enddo

               ! Now compute limited fourth order slopes
               do j = ilo2-1, ihi2+1
                  dq1 = four3rd*dcen(i,j) - sixth*( df(i,j+1) + df(i,j-1) )
                  dqy(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
               enddo
            enddo

            ! Compute slopes in third coordinate direction
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1

                  ! Compute Fromm slope on slab below
                  k = k3d-1
                  dm = two*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = two*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = fourth*(dm+dp)
                  ds = sign( one, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. zero) then
                     dl = sl
                  else
                     dl = zero
                  endif
                  dfm = ds*min(dl,abs(dc))

                  ! Compute Fromm slope on slab above
                  k = k3d+1
                  dm = two*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = two*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = fourth*(dm+dp)
                  ds = sign( one, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. zero) then
                     dl = sl
                  else
                     dl = zero
                  endif
                  dfp = ds*min(dl,abs(dc))

                  ! Compute Fromm slope on current slab
                  k = k3d
                  dm = two*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = two*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = fourth*(dm+dp)
                  ds = sign( one, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. zero) then
                     dl = sl
                  else
                     dl = zero
                  endif

                  ! Now compute limited fourth order slopes
                  dq1 = four3rd*dc - sixth*( dfp + dfm )
                  dqz(i,j,kc,n) = flatn(i,j,k3d)*ds*min(dl,abs(dq1))
               enddo
            enddo
         enddo

      endif

      deallocate(dsgn,dlim,df,dcen)

      end subroutine uslope

end module slope_module
````

