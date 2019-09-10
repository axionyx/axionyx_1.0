
      subroutine ca_axphase(phase,phase_l1,phase_l2,phase_l3,phase_h1,phase_h2,phase_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

      ! use fdm_params_module, only : mindens, meandens

      implicit none

      integer          lo(3), hi(3)
      integer          phase_l1,phase_l2,phase_l3,phase_h1,phase_h2,phase_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision phase(phase_l1:phase_h1,phase_l2:phase_h2,phase_l3:phase_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i, j, k !, n

!       n = 0

!       do k = lo(3), hi(3)                                                                                                                                                          
!          do j = lo(2), hi(2)                                                                                                                                                       
!             do i = lo(1), hi(1)
!                if (dat(i,j,k,1) .gt. (0.75d0*meandens) ) then
!                   n = n+8  !Since the boson star is equally located inside 8 patches
!                endif
!             enddo
!          enddo
!       enddo

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)
!                if (dat(i,j,k,1) .gt. (0.75d0*meandens) ) then
               phase(i,j,k,1)= datan2(dat(i,j,k,3),dat(i,j,k,2)) !/(n*delta(1)*delta(2)*delta(3))
!                else
!                   phase(i,j,k,1)=0.0
!                endif
            enddo
         enddo
      enddo
end subroutine ca_axphase



!-----------------------------------------------------------------------

      subroutine ca_axepot(epot,epot_l1,epot_l2,epot_l3,epot_h1,epot_h2,epot_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

      implicit none

      integer          lo(3), hi(3)
      integer          epot_l1,epot_l2,epot_l3,epot_h1,epot_h2,epot_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision epot(epot_l1:epot_h1,epot_l2:epot_h2,epot_l3:epot_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i, j, k

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)
                  epot(i,j,k,1)= -dat(i,j,k,1)*dat(i,j,k,2)/2
            enddo
         enddo
      enddo

      end subroutine ca_axepot

!-----------------------------------------------------------------------

      subroutine ca_axekin(ekin,ekin_l1,ekin_l2,ekin_l3,ekin_h1,ekin_h2,ekin_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

      use fdm_params_module, only : m_tt,meandens,hbaroverm

      implicit none

      integer          lo(3), hi(3)
      integer          ekin_l1,ekin_l2,ekin_l3,ekin_h1,ekin_h2,ekin_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision ekin(ekin_l1:ekin_h1,ekin_l2:ekin_h2,ekin_l3:ekin_h3,nk)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i, j, k

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)
               ekin(i,j,k,1)= hbaroverm**2 / 2 * &
                              (((-dat(i+2,j,k,2)+8.0d0*dat(i+1,j,k,2)-8.0d0*dat(i-1,j,k,2)+dat(i-2,j,k,2))/(12.0d0*delta(1)))**2 &
                              +((-dat(i+2,j,k,3)+8.0d0*dat(i+1,j,k,3)-8.0d0*dat(i-1,j,k,3)+dat(i-2,j,k,3))/(12.0d0*delta(1)))**2 &
                              +((-dat(i,j+2,k,2)+8.0d0*dat(i,j+1,k,2)-8.0d0*dat(i,j-1,k,2)+dat(i,j-2,k,2))/(12.0d0*delta(2)))**2 &
                              +((-dat(i,j+2,k,3)+8.0d0*dat(i,j+1,k,3)-8.0d0*dat(i,j-1,k,3)+dat(i,j-2,k,3))/(12.0d0*delta(2)))**2 &
                              +((-dat(i,j,k+2,2)+8.0d0*dat(i,j,k+1,2)-8.0d0*dat(i,j,k-1,2)+dat(i,j,k-2,2))/(12.0d0*delta(3)))**2 &
                              +((-dat(i,j,k+2,3)+8.0d0*dat(i,j,k+1,3)-8.0d0*dat(i,j,k-1,3)+dat(i,j,k-2,3))/(12.0d0*delta(3)))**2 )

!                               +((dat(i,j+1,k,2)-dat(i,j-1,k,2))/(2.0d0*delta(2)))**2+((dat(i,j+1,k,3)-dat(i,j-1,k,3))/(2.0d0*delta(2)))**2 &
!                               +((dat(i,j,k+1,2)-dat(i,j,k-1,2))/(2.0d0*delta(3)))**2+((dat(i,j,k+1,3)-dat(i,j,k-1,3))/(2.0d0*delta(3)))**2 )




!                               (((dat(i+1,j,k,2)-dat(i-1,j,k,2))/(2.0d0*delta(1)))**2+((dat(i+1,j,k,3)-dat(i-1,j,k,3))/(2.0d0*delta(1)))**2 &
!                               +((dat(i,j+1,k,2)-dat(i,j-1,k,2))/(2.0d0*delta(2)))**2+((dat(i,j+1,k,3)-dat(i,j-1,k,3))/(2.0d0*delta(2)))**2 &
!                               +((dat(i,j,k+1,2)-dat(i,j,k-1,2))/(2.0d0*delta(3)))**2+((dat(i,j,k+1,3)-dat(i,j,k-1,3))/(2.0d0*delta(3)))**2 )


!                                  (( (-dsqrt(dat(i+2,j,k,1))+8.0d0*dsqrt(dat(i+1,j,k,1))-8.0d0*dsqrt(dat(i-1,j,k,1))+dsqrt(dat(i-2,j,k,1))) / (12.0d0*delta(1)) )**2 &
!                                  +( (-dsqrt(dat(i,j+2,k,1))+8.0d0*dsqrt(dat(i,j+1,k,1))-8.0d0*dsqrt(dat(i,j-1,k,1))+dsqrt(dat(i,j-2,k,1))) / (12.0d0*delta(2)) )**2 &
!                                  +( (-dsqrt(dat(i,j,k+2,1))+8.0d0*dsqrt(dat(i,j,k+1,1))-8.0d0*dsqrt(dat(i,j,k-1,1))+dsqrt(dat(i,j,k-2,1))) / (12.0d0*delta(3)) )**2 )

            enddo
         enddo
      enddo

      end subroutine ca_axekin

!-----------------------------------------------------------------------

      subroutine ca_axekinrho(ekinrho,ekinrho_l1,ekinrho_l2,ekinrho_l3,ekinrho_h1,ekinrho_h2,ekinrho_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

      use fdm_params_module, only : m_tt,hbaroverm,a

      implicit none

      integer          lo(3), hi(3)
      integer          ekinrho_l1,ekinrho_l2,ekinrho_l3,ekinrho_h1,ekinrho_h2,ekinrho_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision ekinrho(ekinrho_l1:ekinrho_h1,ekinrho_l2:ekinrho_h2,ekinrho_l3:ekinrho_h3,nk)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i, j, k

      !Unfortunately, the interpolation argorithm sometimes yields a small negative density
      !if the density is close to zero in a specific area. We can't allow that and set it to zero. 

      do k = lo(3)-2, hi(3)+2                                                                                                                                                          
         do j = lo(2)-2, hi(2)+2                                                                                                                                                       
            do i = lo(1)-2, hi(1)+2
               if (dat(i,j,k,1) .lt. 0.0d0) then
                  dat(i,j,k,1) = 0.0d0
               endif
            enddo
         enddo
      enddo

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)
               ekinrho(i,j,k,1)= hbaroverm**2 / 2 * &
                                 (( (-dsqrt(dat(i+2,j,k,1))+8.0d0*dsqrt(dat(i+1,j,k,1))-8.0d0*dsqrt(dat(i-1,j,k,1))+dsqrt(dat(i-2,j,k,1))) / (12.0d0*delta(1)*a) )**2 &
                                 +( (-dsqrt(dat(i,j+2,k,1))+8.0d0*dsqrt(dat(i,j+1,k,1))-8.0d0*dsqrt(dat(i,j-1,k,1))+dsqrt(dat(i,j-2,k,1))) / (12.0d0*delta(2)*a) )**2 &
                                 +( (-dsqrt(dat(i,j,k+2,1))+8.0d0*dsqrt(dat(i,j,k+1,1))-8.0d0*dsqrt(dat(i,j,k-1,1))+dsqrt(dat(i,j,k-2,1))) / (12.0d0*delta(3)*a) )**2 )
            enddo
         enddo
      enddo

      end subroutine ca_axekinrho

!-----------------------------------------------------------------------

      subroutine ca_axekinv(ekinv,ekinv_l1,ekinv_l2,ekinv_l3,ekinv_h1,ekinv_h2,ekinv_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

      use fdm_params_module, only : m_tt, hbaroverm,a!, mindens
      use fundamental_constants_module

      implicit none

      integer          lo(3), hi(3)
      integer          ekinv_l1,ekinv_l2,ekinv_l3,ekinv_h1,ekinv_h2,ekinv_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt, diff(3)
      double precision ekinv(ekinv_l1:ekinv_h1,ekinv_l2:ekinv_h2,ekinv_l3:ekinv_h3,nk)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      !double precision phase(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3)
      double precision phase(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
      integer    level, grid_no

      integer i, j, k

      if (nc .eq. 3) then !We solve AxDens,AxIm,AxRe 

      do k = lo(3)-1, hi(3)+1                                                                                                                                                          
         do j = lo(2)-1, hi(2)+1                                                                                                                                                       
            do i = lo(1)-1, hi(1)+1
               if (dat(i,j,k,1) .ne. 0.0d0) then
                  phase(i,j,k) = datan2(dat(i,j,k,3),dat(i,j,k,2))
               else
                  phase(i,j,k)=0.0
               endif
            enddo
         enddo
      enddo

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)

               !Have to consider the cases where the phase jumps from close to -pi to pi and vice versa
               diff(1) = min( dabs(phase(i+1,j,k)-phase(i-1,j,k)) , dabs(phase(i+1,j,k)-phase(i-1,j,k)-2.0*PI) , dabs(phase(i+1,j,k)-phase(i-1,j,k)+2.0*PI) ) 
               diff(2) = min( dabs(phase(i,j+1,k)-phase(i,j-1,k)) , dabs(phase(i,j+1,k)-phase(i,j-1,k)-2.0*PI) , dabs(phase(i,j+1,k)-phase(i,j-1,k)+2.0*PI) ) 
               diff(3) = min( dabs(phase(i,j,k+1)-phase(i,j,k-1)) , dabs(phase(i,j,k+1)-phase(i,j,k-1)-2.0*PI) , dabs(phase(i,j,k+1)-phase(i,j,k-1)+2.0*PI) ) 

               ekinv(i,j,k,1)= hbaroverm**2 / 2 * dat(i,j,k,1) * &
                               (( diff(1) / (2*delta(1)*a) )**2 &
                               +( diff(2) / (2*delta(2)*a) )**2 &
                               +( diff(3) / (2*delta(3)*a) )**2)

               ! !Have to consider the cases where the phase jumps from close to -pi to pi and vice versa
               ! diff(1) = min( dabs(phase(i+1,j,k)-phase(i,j,k)) , dabs(phase(i+1,j,k)-phase(i,j,k)-2.0*PI) , dabs(phase(i+1,j,k)-phase(i,j,k)+2.0*PI) ) 
               ! diff(2) = min( dabs(phase(i,j+1,k)-phase(i,j,k)) , dabs(phase(i,j+1,k)-phase(i,j,k)-2.0*PI) , dabs(phase(i,j+1,k)-phase(i,j,k)+2.0*PI) ) 
               ! diff(3) = min( dabs(phase(i,j,k+1)-phase(i,j,k)) , dabs(phase(i,j,k+1)-phase(i,j,k)-2.0*PI) , dabs(phase(i,j,k+1)-phase(i,j,k)+2.0*PI) ) 

               ! ekinv(i,j,k,1)= hbaroverm**2 / 2 * ( dat(i,j,k,1)+dat(i+1,j,k,1)+dat(i,j+1,k,1)+dat(i,j,k+1,1)+dat(i+1,j+1,k,1)+dat(i,j+1,k+1,1)+dat(i+1,j,k+1,1)+dat(i+1,j+1,k+1,1) )/8 * &
               !                 (( diff(1) / delta(1) )**2 &
               !                 +( diff(2) / delta(2) )**2 &
               !                 +( diff(3) / delta(3) )**2)
            enddo
         enddo
      enddo

      else !We solve UAXDENS, UAXMOMX, UAXMOMY, UAXMOMZ
         
      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)

               if (dat(i,j,k,1) .eq. 0.d0 ) then
                  ekinv(i,j,k,1) = 0.d0
               else
                  ekinv(i,j,k,1)= ( dat(i,j,k,2)**2 + dat(i,j,k,3)**2 + dat(i,j,k,4)**2 )/dat(i,j,k,1)/2.d0
               endif

            enddo
         enddo
      enddo
                  
      endif

      end subroutine ca_axekinv


!-----------------------------------------------------------------------

      subroutine ca_axvel(ekinv,ekinv_l1,ekinv_l2,ekinv_l3,ekinv_h1,ekinv_h2,ekinv_h3,nk, &
                          dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                          lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

      use fdm_params_module, only : m_tt, hbaroverm!, mindens
      use fundamental_constants_module

      implicit none

      integer          lo(3), hi(3)
      integer          ekinv_l1,ekinv_l2,ekinv_l3,ekinv_h1,ekinv_h2,ekinv_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt, diff(3)
      double precision ekinv(ekinv_l1:ekinv_h1,ekinv_l2:ekinv_h2,ekinv_l3:ekinv_h3,nk)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      double precision phase(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3)
      integer    level, grid_no

      integer i, j, k

      do k = lo(3)-1, hi(3)+1                                                                                                                                                          
         do j = lo(2)-1, hi(2)+1                                                                                                                                                       
            do i = lo(1)-1, hi(1)+1
               if (dat(i,j,k,1) .ne. 0.0d0) then
                  phase(i,j,k) = datan2(dat(i,j,k,3),dat(i,j,k,2))
               else
                  phase(i,j,k)=0.0
               endif
            enddo
         enddo
      enddo

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)

               !Have to consider the cases where the phase jumps from close to -pi to pi and vice versa
               diff(1) = min( dabs(phase(i+1,j,k)-phase(i-1,j,k)) , dabs(phase(i+1,j,k)-phase(i-1,j,k)-2.0*PI) , dabs(phase(i+1,j,k)-phase(i-1,j,k)+2.0*PI) ) 
               diff(2) = min( dabs(phase(i,j+1,k)-phase(i,j-1,k)) , dabs(phase(i,j+1,k)-phase(i,j-1,k)-2.0*PI) , dabs(phase(i,j+1,k)-phase(i,j-1,k)+2.0*PI) ) 
               diff(3) = min( dabs(phase(i,j,k+1)-phase(i,j,k-1)) , dabs(phase(i,j,k+1)-phase(i,j,k-1)-2.0*PI) , dabs(phase(i,j,k+1)-phase(i,j,k-1)+2.0*PI) ) 

               ekinv(i,j,k,1)= hbaroverm * sqrt(&
                                ( diff(1) / (2*delta(1)) )**2 &
                               +( diff(2) / (2*delta(2)) )**2 &
                               +( diff(3) / (2*delta(3)) )**2)

            enddo
         enddo
      enddo

      end subroutine ca_axvel


 !-----------------------------------------------------------------------
 
      subroutine ca_axangmom_x(angmom_x,angmom_x_l1,angmom_x_l2,angmom_x_l3,angmom_x_h1,angmom_x_h2,angmom_x_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

      use fdm_params_module, only : m_tt,meandens,hbaroverm
      use fundamental_constants_module

      implicit none

      integer          lo(3), hi(3)
      integer          angmom_x_l1,angmom_x_l2,angmom_x_l3,angmom_x_h1,angmom_x_h2,angmom_x_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision angmom_x(angmom_x_l1:angmom_x_h1,angmom_x_l2:angmom_x_h2,angmom_x_l3:angmom_x_h3,nk)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      ! double precision phase(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3)
      integer    level, grid_no

      integer i, j, k

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)

               angmom_x(i,j,k,1) = hbaroverm*meandens*( -dat(i,j,k,2)*(dat(i,j,k+1,3)-dat(i,j,k-1,3))/(2.0d0*delta(3))*(j*delta(2)+xlo(2)) &
                                                        +dat(i,j,k,3)*(dat(i,j,k+1,2)-dat(i,j,k-1,2))/(2.0d0*delta(3))*(j*delta(2)+xlo(2)) &   
                                                        +dat(i,j,k,2)*(dat(i,j+1,k,3)-dat(i,j-1,k,3))/(2.0d0*delta(2))*(k*delta(3)+xlo(3)) &  
                                                        -dat(i,j,k,3)*(dat(i,j+1,k,2)-dat(i,j-1,k,2))/(2.0d0*delta(2))*(k*delta(3)+xlo(3)) )

            enddo 
         enddo
      enddo

      end subroutine ca_axangmom_x

!-----------------------------------------------------------------------

      subroutine ca_axangmom_y(angmom_y,angmom_y_l1,angmom_y_l2,angmom_y_l3,angmom_y_h1,angmom_y_h2,angmom_y_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

      use fdm_params_module, only : m_tt,meandens,hbaroverm
      use fundamental_constants_module

      implicit none

      integer          lo(3), hi(3)
      integer          angmom_y_l1,angmom_y_l2,angmom_y_l3,angmom_y_h1,angmom_y_h2,angmom_y_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision angmom_y(angmom_y_l1:angmom_y_h1,angmom_y_l2:angmom_y_h2,angmom_y_l3:angmom_y_h3,nk)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      ! double precision phase(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3)
      integer    level, grid_no

      integer i, j, k

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)

               angmom_y(i,j,k,1) = hbaroverm*meandens*( -dat(i,j,k,2)*(dat(i+1,j,k,3)-dat(i-1,j,k,3))/(2.0d0*delta(1))*(k*delta(3)+xlo(3)) &
                                                        +dat(i,j,k,3)*(dat(i+1,j,k,2)-dat(i-1,j,k,2))/(2.0d0*delta(1))*(k*delta(3)+xlo(3)) &   
                                                        +dat(i,j,k,2)*(dat(i,j,k+1,3)-dat(i,j,k-1,3))/(2.0d0*delta(3))*(i*delta(1)+xlo(1)) &  
                                                        -dat(i,j,k,3)*(dat(i,j,k+1,2)-dat(i,j,k-1,2))/(2.0d0*delta(3))*(i*delta(1)+xlo(1)) )

            enddo
         enddo
      enddo

      end subroutine ca_axangmom_y

!-----------------------------------------------------------------------

      subroutine ca_axangmom_z(angmom_z,angmom_z_l1,angmom_z_l2,angmom_z_l3,angmom_z_h1,angmom_z_h2,angmom_z_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

      use fdm_params_module, only : m_tt,meandens,hbaroverm
      use fundamental_constants_module

      implicit none

      integer          lo(3), hi(3)
      integer          angmom_z_l1,angmom_z_l2,angmom_z_l3,angmom_z_h1,angmom_z_h2,angmom_z_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision angmom_z(angmom_z_l1:angmom_z_h1,angmom_z_l2:angmom_z_h2,angmom_z_l3:angmom_z_h3,nk)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      ! double precision phase(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3)
      integer    level, grid_no

      integer i, j, k

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)

               angmom_z(i,j,k,1) = hbaroverm*meandens*( -dat(i,j,k,2)*(dat(i,j+1,k,3)-dat(i,j-1,k,3))/(2.0d0*delta(2))*(i*delta(1)+xlo(1)) &
                                                        +dat(i,j,k,3)*(dat(i,j+1,k,2)-dat(i,j-1,k,2))/(2.0d0*delta(2))*(i*delta(1)+xlo(1)) &   
                                                        +dat(i,j,k,2)*(dat(i+1,j,k,3)-dat(i-1,j,k,3))/(2.0d0*delta(1))*(j*delta(2)+xlo(2)) &  
                                                        -dat(i,j,k,3)*(dat(i+1,j,k,2)-dat(i-1,j,k,2))/(2.0d0*delta(1))*(j*delta(2)+xlo(2)) )

            enddo
         enddo
      enddo

      end subroutine ca_axangmom_z

!-----------------------------------------------------------------------


      subroutine ca_lohnererror(err,err_l1,err_l2,err_l3,err_h1,err_h2,err_h3,nk, &
                                dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                                lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine calculates the error estimator for the fdm velocity (AxRe,AxIm) 
      ! according to R. Loehner (1987) eq.4
      !

      use fdm_params_module, only : epsilon_L, mindens

      implicit none

      integer          lo(3), hi(3)
      integer          err_l1,err_l2,err_l3,err_h1,err_h2,err_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision err(err_l1:err_h1,err_l2:err_h2,err_l3:err_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      double precision top, bottom, eps
      integer i, j, k

      eps = 0.2

      do k = lo(3), hi(3)                                                                                                                                      
         do j = lo(2), hi(2)                                                                                                                                              
            do i = lo(1), hi(1)
               
               err(i,j,k,1) = 0.0
              
               ! print *, "dat",dat(i,j,k,1),dat(i,j,k,2)

               if ( (dat(i,j,k,1).ge.mindens) .and. (dat(i,j,k,2).ne.0.0) ) then

               top =       (dat(i+1,j  ,k  ,2)-dat(i  ,j  ,k  ,2)-dat(i  ,j  ,k  ,2)+dat(i-1,j  ,k  ,2))**2 + &
                           (dat(i  ,j+1,k  ,2)-dat(i  ,j  ,k  ,2)-dat(i  ,j  ,k  ,2)+dat(i  ,j-1,k  ,2))**2 + &
                           (dat(i  ,j  ,k+1,2)-dat(i  ,j  ,k  ,2)-dat(i  ,j  ,k  ,2)+dat(i  ,j  ,k-1,2))**2 + &
                     0.125*(dat(i+1,j+1,k  ,2)-dat(i+1,j-1,k  ,2)-dat(i-1,j+1,k  ,2)+dat(i-1,j-1,k  ,2))**2 + &
                     0.125*(dat(i+1,j  ,k+1,2)-dat(i+1,j  ,k-1,2)-dat(i-1,j  ,k+1,2)+dat(i-1,j  ,k-1,2))**2 + &
                     0.125*(dat(i  ,j+1,k+1,2)-dat(i  ,j+1,k-1,2)-dat(i  ,j-1,k+1,2)+dat(i  ,j-1,k-1,2))**2

               bottom =      (abs(dat(i+1,j  ,k  ,2) -    dat(i  ,j  ,k  ,2))+abs(dat(i  ,j  ,k  ,2) -    dat(i-1,j  ,k  ,2))      + &
                         eps*(abs(dat(i+1,j  ,k  ,2))+abs(dat(i  ,j  ,k  ,2))+abs(dat(i  ,j  ,k  ,2))+abs(dat(i-1,j  ,k  ,2))))**2 + &
                             (abs(dat(i  ,j+1,k  ,2) -    dat(i  ,j  ,k  ,2))+abs(dat(i  ,j  ,k  ,2) -    dat(i  ,j-1,k  ,2))      + &
                         eps*(abs(dat(i  ,j+1,k  ,2))+abs(dat(i  ,j  ,k  ,2))+abs(dat(i  ,j  ,k  ,2))+abs(dat(i  ,j-1,k  ,2))))**2 + &
                             (abs(dat(i  ,j  ,k+1,2) -    dat(i  ,j  ,k  ,2))+abs(dat(i  ,j  ,k  ,2) -    dat(i  ,j  ,k-1,2))      + &
                         eps*(abs(dat(i  ,j  ,k+1,2))+abs(dat(i  ,j  ,k  ,2))+abs(dat(i  ,j  ,k  ,2))+abs(dat(i  ,j  ,k-1,2))))**2 + &
                        (0.5*(abs(dat(i+1,j+1,k  ,2) -    dat(i-1,j+1,k  ,2))+abs(dat(i+1,j-1,k  ,2) -    dat(i-1,j-1,k  ,2)))     + &
                         eps*(abs(dat(i+1,j+1,k  ,2))+abs(dat(i-1,j+1,k  ,2))+abs(dat(i+1,j-1,k  ,2))+abs(dat(i-1,j-1,k  ,2))))**2 + &
                        (0.5*(abs(dat(i+1,j+1,k  ,2) -    dat(i+1,j-1,k  ,2))+abs(dat(i-1,j+1,k  ,2) -    dat(i-1,j-1,k  ,2)))     + &
                         eps*(abs(dat(i+1,j+1,k  ,2))+abs(dat(i+1,j-1,k  ,2))+abs(dat(i-1,j+1,k  ,2))+abs(dat(i-1,j-1,k  ,2))))**2 + &
                        (0.5*(abs(dat(i+1,j  ,k+1,2) -    dat(i-1,j  ,k+1,2))+abs(dat(i+1,j  ,k-1,2) -    dat(i-1,j  ,k-1,2)))     + &
                         eps*(abs(dat(i+1,j  ,k+1,2))+abs(dat(i-1,j  ,k+1,2))+abs(dat(i+1,j  ,k-1,2))+abs(dat(i-1,j  ,k-1,2))))**2 + &
                        (0.5*(abs(dat(i+1,j  ,k+1,2) -    dat(i+1,j  ,k-1,2))+abs(dat(i-1,j  ,k+1,2) -    dat(i-1,j  ,k-1,2)))     + &
                         eps*(abs(dat(i+1,j  ,k+1,2))+abs(dat(i+1,j  ,k-1,2))+abs(dat(i-1,j  ,k+1,2))+abs(dat(i-1,j  ,k-1,2))))**2 + &
                        (0.5*(abs(dat(i  ,j+1,k+1,2) -    dat(i  ,j-1,k+1,2))+abs(dat(i  ,j+1,k-1,2) -    dat(i  ,j-1,k-1,2)))     + &
                         eps*(abs(dat(i  ,j+1,k+1,2))+abs(dat(i  ,j-1,k+1,2))+abs(dat(i  ,j+1,k-1,2))+abs(dat(i  ,j-1,k-1,2))))**2 + &
                        (0.5*(abs(dat(i  ,j+1,k+1,2) -    dat(i  ,j+1,k-1,2))+abs(dat(i  ,j-1,k+1,2) -    dat(i  ,j-1,k-1,2)))     + &
                         eps*(abs(dat(i  ,j+1,k+1,2))+abs(dat(i  ,j+1,k-1,2))+abs(dat(i  ,j-1,k+1,2))+abs(dat(i  ,j-1,k-1,2))))**2

               err(i,j,k,1) = sqrt(top/bottom)
               endif

               if ( (dat(i,j,k,1).ge.mindens) .and. (dat(i,j,k,3).ne.0.0) ) then

               top =       (dat(i+1,j  ,k  ,3)-dat(i  ,j  ,k  ,3)-dat(i  ,j  ,k  ,3)+dat(i-1,j  ,k  ,3))**2 + &
                           (dat(i  ,j+1,k  ,3)-dat(i  ,j  ,k  ,3)-dat(i  ,j  ,k  ,3)+dat(i  ,j-1,k  ,3))**2 + &
                           (dat(i  ,j  ,k+1,3)-dat(i  ,j  ,k  ,3)-dat(i  ,j  ,k  ,3)+dat(i  ,j  ,k-1,3))**2 + &
                     0.125*(dat(i+1,j+1,k  ,3)-dat(i+1,j-1,k  ,3)-dat(i-1,j+1,k  ,3)+dat(i-1,j-1,k  ,3))**2 + &
                     0.125*(dat(i+1,j  ,k+1,3)-dat(i+1,j  ,k-1,3)-dat(i-1,j  ,k+1,3)+dat(i-1,j  ,k-1,3))**2 + &
                     0.125*(dat(i  ,j+1,k+1,3)-dat(i  ,j+1,k-1,3)-dat(i  ,j-1,k+1,3)+dat(i  ,j-1,k-1,3))**2

               bottom =      (abs(dat(i+1,j  ,k  ,3) -    dat(i  ,j  ,k  ,3))+abs(dat(i  ,j  ,k  ,3) -    dat(i-1,j  ,k  ,3))      + &
                         eps*(abs(dat(i+1,j  ,k  ,3))+abs(dat(i  ,j  ,k  ,3))+abs(dat(i  ,j  ,k  ,3))+abs(dat(i-1,j  ,k  ,3))))**2 + &
                             (abs(dat(i  ,j+1,k  ,3) -    dat(i  ,j  ,k  ,3))+abs(dat(i  ,j  ,k  ,3) -    dat(i  ,j-1,k  ,3))      + &
                         eps*(abs(dat(i  ,j+1,k  ,3))+abs(dat(i  ,j  ,k  ,3))+abs(dat(i  ,j  ,k  ,3))+abs(dat(i  ,j-1,k  ,3))))**2 + &
                             (abs(dat(i  ,j  ,k+1,3) -    dat(i  ,j  ,k  ,3))+abs(dat(i  ,j  ,k  ,3) -    dat(i  ,j  ,k-1,3))      + &
                         eps*(abs(dat(i  ,j  ,k+1,3))+abs(dat(i  ,j  ,k  ,3))+abs(dat(i  ,j  ,k  ,3))+abs(dat(i  ,j  ,k-1,3))))**2 + &
                        (0.5*(abs(dat(i+1,j+1,k  ,3) -    dat(i-1,j+1,k  ,3))+abs(dat(i+1,j-1,k  ,3) -    dat(i-1,j-1,k  ,3)))     + &
                         eps*(abs(dat(i+1,j+1,k  ,3))+abs(dat(i-1,j+1,k  ,3))+abs(dat(i+1,j-1,k  ,3))+abs(dat(i-1,j-1,k  ,3))))**2 + &
                        (0.5*(abs(dat(i+1,j+1,k  ,3) -    dat(i+1,j-1,k  ,3))+abs(dat(i-1,j+1,k  ,3) -    dat(i-1,j-1,k  ,3)))     + &
                         eps*(abs(dat(i+1,j+1,k  ,3))+abs(dat(i+1,j-1,k  ,3))+abs(dat(i-1,j+1,k  ,3))+abs(dat(i-1,j-1,k  ,3))))**2 + &
                        (0.5*(abs(dat(i+1,j  ,k+1,3) -    dat(i-1,j  ,k+1,3))+abs(dat(i+1,j  ,k-1,3) -    dat(i-1,j  ,k-1,3)))     + &
                         eps*(abs(dat(i+1,j  ,k+1,3))+abs(dat(i-1,j  ,k+1,3))+abs(dat(i+1,j  ,k-1,3))+abs(dat(i-1,j  ,k-1,3))))**2 + &
                        (0.5*(abs(dat(i+1,j  ,k+1,3) -    dat(i+1,j  ,k-1,3))+abs(dat(i-1,j  ,k+1,3) -    dat(i-1,j  ,k-1,3)))     + &
                         eps*(abs(dat(i+1,j  ,k+1,3))+abs(dat(i+1,j  ,k-1,3))+abs(dat(i-1,j  ,k+1,3))+abs(dat(i-1,j  ,k-1,3))))**2 + &
                        (0.5*(abs(dat(i  ,j+1,k+1,3) -    dat(i  ,j-1,k+1,3))+abs(dat(i  ,j+1,k-1,3) -    dat(i  ,j-1,k-1,3)))     + &
                         eps*(abs(dat(i  ,j+1,k+1,3))+abs(dat(i  ,j-1,k+1,3))+abs(dat(i  ,j+1,k-1,3))+abs(dat(i  ,j-1,k-1,3))))**2 + &
                        (0.5*(abs(dat(i  ,j+1,k+1,3) -    dat(i  ,j+1,k-1,3))+abs(dat(i  ,j-1,k+1,3) -    dat(i  ,j-1,k-1,3)))     + &
                         eps*(abs(dat(i  ,j+1,k+1,3))+abs(dat(i  ,j+1,k-1,3))+abs(dat(i  ,j-1,k+1,3))+abs(dat(i  ,j-1,k-1,3))))**2

                  err(i,j,k,1) = max(err(i,j,k,1),sqrt(top/bottom))

                  endif

                  ! print *, 'err',err(i,j,k,1)

            enddo
         enddo
      enddo

      end subroutine ca_lohnererror

!-----------------------------------------------------------------------


      subroutine ca_dererrx(err,err_l1,err_l2,err_l3,err_h1,err_h2,err_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine calculates the error estimator for the fdm velocity (AxRe,AxIm) 
      ! according to R. Loehner (1987) eq.4
      !

      use fdm_params_module, only : epsilon_L, mindens

      implicit none

      integer          lo(3), hi(3)
      integer          err_l1,err_l2,err_l3,err_h1,err_h2,err_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision err(err_l1:err_h1,err_l2:err_h2,err_l3:err_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      double precision secondder, firstder, meanvalue
      integer i, j, k

      do k = lo(3), hi(3)                                                                                                                                      
         do j = lo(2), hi(2)                                                                                                                                              
            do i = lo(1), hi(1)

               secondder = dabs(dat(i+1,j,k,1)-2.*dat(i,j,k,1)+dat(i-1,j,k,1))
               firstder  = dabs(dat(i+1,j,k,1)-dat(i,j,k,1)) +dabs(dat(i,j,k,1)-dat(i-1,j,k,1))
               meanvalue = dabs(dat(i+1,j,k,1)+2.*dat(i,j,k,1)+dat(i-1,j,k,1))                                                 
               if ( ( firstder .gt. 1.0d-100 ) .and. ( dabs(dat(i,j,k,2)) .gt. mindens ) ) then
                  err(i,j,k,1) = secondder/(firstder+epsilon_L*meanvalue)
               else
                  err(i,j,k,1) = 0.0d0
               endif

            enddo
         enddo
      enddo

      end subroutine ca_dererrx

!-----------------------------------------------------------------------

      subroutine ca_dererry(err,err_l1,err_l2,err_l3,err_h1,err_h2,err_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine calculates the error estimator for the fdm velocity (AxRe,AxIm) 
      ! according to R. Loehner (1987) eq.4
      !

      use fdm_params_module, only : epsilon_L, mindens

      implicit none

      integer          lo(3), hi(3)
      integer          err_l1,err_l2,err_l3,err_h1,err_h2,err_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision err(err_l1:err_h1,err_l2:err_h2,err_l3:err_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      double precision secondder, firstder, meanvalue
      integer i, j, k

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)

               secondder = dabs(dat(i,j+1,k,1)-2.*dat(i,j,k,1)+dat(i,j-1,k,1))
               firstder  = dabs(dat(i,j+1,k,1)-dat(i,j,k,1)) +dabs(dat(i,j,k,1)-dat(i,j-1,k,1))
               meanvalue = dabs(dat(i,j+1,k,1)+2.*dat(i,j,k,1)+dat(i,j-1,k,1))
               if ( ( firstder .gt. 1.0d-100 ) .and. ( dabs(dat(i,j,k,2)) .gt. mindens ) ) then
                  err(i,j,k,1) = secondder/(firstder+epsilon_L*meanvalue)
               else
                  err(i,j,k,1) = 0.0d0
               endif

            enddo
         enddo
      enddo

      end subroutine ca_dererry

!-----------------------------------------------------------------------

      subroutine ca_dererrz(err,err_l1,err_l2,err_l3,err_h1,err_h2,err_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine calculates the error estimator for the fdm velocity (AxRe,AxIm) 
      ! according to R. Loehner (1987) eq.4
      !

      use fdm_params_module, only : epsilon_L, mindens

      implicit none

      integer          lo(3), hi(3)
      integer          err_l1,err_l2,err_l3,err_h1,err_h2,err_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision err(err_l1:err_h1,err_l2:err_h2,err_l3:err_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      double precision secondder, firstder, meanvalue
      integer i, j, k

      do k = lo(3), hi(3)                                                                                                                                                          
         do j = lo(2), hi(2)                                                                                                                                                       
            do i = lo(1), hi(1)

               secondder = dabs(dat(i,j,k+1,1)-2.*dat(i,j,k,1)+dat(i,j,k-1,1))
               firstder  = dabs(dat(i,j,k+1,1)-dat(i,j,k,1)) +dabs(dat(i,j,k,1)-dat(i,j,k-1,1))
               meanvalue = dabs(dat(i,j,k+1,1)+2.*dat(i,j,k,1)+dat(i,j,k-1,1))                                                 
               if ( ( firstder .gt. 1.0d-100 ) .and. ( dabs(dat(i,j,k,2)) .gt. mindens ) ) then
                  err(i,j,k,1) = secondder/(firstder+epsilon_L*meanvalue)
               else
                  err(i,j,k,1) = 0.0d0
               endif

            enddo
         enddo
      enddo

      end subroutine ca_dererrz
