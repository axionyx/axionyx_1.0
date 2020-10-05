
! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on overdensity
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nc        => number of components in density array
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_overdensity(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                 set,clear, &
                                 den,denl1,denl2,denl3,denh1,denh2,denh3, &
                                 lo,hi,nc,domlo,domhi,delta,level,avg_den)

      use probdata_module
      use comoving_module, only : comoving_h, comoving_OmAx
      implicit none

      integer set, clear, nc, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer denl1,denl2,denl3,denh1,denh2,denh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision den(denl1:denh1,denl2:denh2,denl3:denh3,nc)
      double precision delta(3), avg_den, critdens

      double precision :: over_den

      integer i, j, k

      over_den = 2.0d0**(3*(level+1)) * avg_den
      ! critdens = 2.775d11 * 0.7d0**2 * comoving_OmAx !3 * (100)**2 *0.7d0**2 *  comoving_OmAx / ( 8 * PI * Gconst)                                                                                                                                                          
      ! over_den = 2.0d0**(3*(level+1))*critdens
!     Tag on regions of overdensity
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( den(i,j,k,1) .gt. over_den ) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
      enddo

      end subroutine tag_overdensity

      subroutine tag_axvel(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              set,clear, &
                              vel,vell1,vell2,vell3,velh1,velh2,velh3, &
                              lo,hi,nc,domlo,domhi,delta,xlo,problo,time,level)

      use probdata_module
      use fdm_params_module, only : critvalue

      implicit none

      integer set, clear, nc, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer vell1,vell2,vell3,velh1,velh2,velh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision vel(vell1:velh1,vell2:velh2,vell3:velh3,nc)
      double precision delta(3), avg_den

      double precision xlo(3), problo(3),time
      integer i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( vel(i,j,k,1) .gt. critvalue ) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
      enddo

      end subroutine tag_axvel

      subroutine tag_center(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                            set,clear, &
                            den,denl1,denl2,denl3,denh1,denh2,denh3, &
                            lo,hi,nc,domlo,domhi,delta,level,avg_den)

      use probdata_module
      use comoving_module, only : comoving_h, comoving_OmAx
      implicit none

      integer set, clear, nc, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer denl1,denl2,denl3,denh1,denh2,denh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision den(denl1:denh1,denl2:denh2,denl3:denh3,nc)
      double precision delta(3), avg_den, critdens

      integer i, j, k
!       double precision low(3), high(3)
      double precision r

!       do i = 1,3
!          low(i)  = (1.d0-0.5d0**(level+1)) * 0.5d0 * domhi(i) !*0.25d0
!          high(i) = (1.d0+0.5d0**(level+1)) * 0.5d0 * domhi(i) !*0.75d0
!       enddo

!       do k = lo(3), hi(3)
!          do j = lo(2), hi(2)
!             do i = lo(1), hi(1)
!                if ( i .lt. high(1) .and. i .gt. low(1) ) then
!                   if ( j .lt. high(2) .and. j .gt. low(2) ) then
!                      if ( k .lt. high(3) .and. k .gt. low(3) ) then
!                         tag(i,j,k) = set
!                      endif
!                   endif
!                endif
!             enddo
!          enddo
!       enddo

      r = (5-level)/16.d0*domhi(1)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( dsqrt( (i-0.5d0*domhi(1))**2 + (j-0.5d0*domhi(2))**2 + (k-0.5d0*domhi(3))**2 ) .lt. r ) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
      enddo

      end subroutine tag_center


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on overdensity
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nc        => number of components in density array
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_lohner(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                            set,clear, &
                            err,errl1,errl2,errl3,errh1,errh2,errh3, &
                            lo,hi,nc,domlo,domhi,delta,level)

      use probdata_module
      implicit none

      integer set, clear, nc, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer errl1,errl2,errl3,errh1,errh2,errh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision err(errl1:errh1,errl2:errh2,errl3:errh3,nc)
      double precision delta(3), avg_den, critdens, maxi

      integer i, j, k

      maxi=0.0
!     Tag on regions where the Lohner criterion is not met.
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( err(i,j,k,1) .gt. 0.1 ) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
      enddo

      end subroutine tag_lohner
