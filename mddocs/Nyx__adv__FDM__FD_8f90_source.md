
# File Nyx\_adv\_FDM\_FD.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Src\_3d**](dir_723248e6e98dc7cb10ec13b7569a328c.md) **>** [**Nyx\_adv\_FDM\_FD.f90**](Nyx__adv__FDM__FD_8f90.md)

[Go to the documentation of this file.](Nyx__adv__FDM__FD_8f90.md) 


````cpp
      subroutine tdma(n,a,b,c,d,x)

!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations

      implicit none

      integer, intent(in) :: n
      double complex, intent(in) ::                  a(n), c(n)
      double complex, intent(inout), dimension(n) :: b, d
      double complex, intent(out) ::                 x(n)
      !  --- Local variables ---
      integer :: i
      double complex :: q

      !  --- Elimination ---
      do i = 2,n
         q = a(i)/b(i - 1)
         b(i) = b(i) - c(i - 1)*q
         d(i) = d(i) - d(i - 1)*q
      enddo

      ! --- Backsubstitution ---
      q = d(n)/b(n)
      x(n) = q
      do i = n - 1,1,-1
         q = (d(i) - c(i)*q)/b(i)
         x(i) = q
      enddo

      return

      end subroutine tdma

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_estdtax(uin,  uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                    uout, uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                    dt_old, dt_new, delta, maxchange, x)

      use meth_params_module, only : uaxre, uaxim
      use fundamental_constants_module
      !use fdm_params_module, only: ax_maxchange,ax_x,ax_y,ax_z

      implicit none

      integer          uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer          uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer          i, j, k

      ! double precision  uin(  uin_l1:uin_h1,   uin_l2:uin_h2,   uin_l3:uin_h3,  NAXVAR)
      ! double precision uout( uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, NAXVAR)
      double precision  uin(  uin_l1:uin_h1,   uin_l2:uin_h2,   uin_l3:uin_h3, 1)
      double precision uout( uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, 1)
      double precision dt_old, dt_new
      double precision ampin, ampout, maxchange_temp, maxchange 
      double precision delta(3), x(3)

      dt_new = dt_old
      maxchange_temp = 0.0d0

      !iterate over physical domain and not over ghost cells
      do k = uout_l3+1, uout_h3-1
        do j = uout_l2+1, uout_h2-1
          do i = uout_l1+1, uout_h1-1

             ampin = sqrt(uin(i,j,k, uaxre)**2 + uin(i,j,k, uaxim)**2)
             ampout = sqrt(uout(i,j,k, uaxre)**2 + uout(i,j,k, uaxim)**2)
             ! ampin = uin(i,j,k,1)
             ! ampout = uout(i,j,k,1)
             ! phasein = datan2(uin(i,j,k, UAXIM), uin(i,j,k, UAXRE))
             ! phaseout = datan2(uout(i,j,k, UAXIM), uout(i,j,k, UAXRE))
             ! if (abs(phasein - phaseout) .gt. 1.9*PI .or. ampin .lt. 1.0d-5) then 
             !    phasein = 1.
             !    phaseout = 1.
             !    !print *, 'oops'
             ! endif
             maxchange_temp = max(maxchange_temp, dabs( ampin/ampout - 1.0 ))!, &
                              !dabs( (phasein - phaseout)/(2.*PI) )) !OLD VERSIONallow phasechange to be 10x larger than ampchange
             
             ! if (maxchange_temp .gt. maxchange) then
             !    maxchange = maxchange_temp
             !    x(1) = i*delta(1)
             !    x(2) = j*delta(2)
             !    x(3) = k*delta(3)
             ! endif

          enddo
        enddo
      enddo
!      print *, maxchange_temp
      if (maxchange_temp .gt. 1.d-1) then
         dt_new = 1.d-1 * dt_old / maxchange_temp
      else 
         dt_new = 2.*dt_old
      endif

      return

      end subroutine fort_estdtax

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_divvel(uin,  uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                             uout, uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                             vdx, nvar, delta, ProbLo,ProbHi &
                             )

      use fdm_params_module, only : hbaroverm

      implicit none

      integer          uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer          uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer          i, j, k
      integer          vdx, nvar

      double precision uin(   uin_l1:uin_h1,  uin_l2:uin_h2,   uin_l3:uin_h3,  nvar)
      double precision uout(uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, 1) 
      double precision delvx, delvy, delvz
      double precision delta(3), ProbLo(3),ProbHi(3),center(3)

      do i = 1,3
         center(i) = (probhi(i)-problo(i))/2.d0
      end do      

      !hbar/m expressed in Nyx units [Mpc km/s]
      ! hbaroverm = 0.01917152d0 / m_tt


      ! do k = uout_l3, uout_h3
      !   do j = uout_l2, uout_h2
      !     do i = uout_l1, uout_h1

      !        r = dsqrt( (i*delta(1)-center(1))**2 + (j*delta(2)-center(2))**2 + (k*delta(3)-center(3))**2 )
      !        if (r .lt. 1.d-6) then 
      !            r = 1.d-6
      !        endif
      !        uout(i,j,k,1) = exp(-r)/r
      !        ! uout(i,j,k,1) = mod(i+j+k,2)*2-1
      !        ! uout(i,j,k,1) = mod(i,2)*2-1

      !     enddo
      !   enddo
      ! enddo

      
      ! need 1 ghost zone -> uout_l = uin_l+1, uout_h = uin_h-1
      do k = uout_l3, uout_h3
        do j = uout_l2, uout_h2
          do i = uout_l1, uout_h1

             delvx = ( uin(i+1,j,k,vdx  ) - uin(i-1,j,k,vdx  ) ) / (2.*delta(1))
             delvy = ( uin(i,j+1,k,vdx+1) - uin(i,j-1,k,vdx+1) ) / (2.*delta(2))
             delvz = ( uin(i,j,k+1,vdx+2) - uin(i,j,k-1,vdx+2) ) / (2.*delta(3))
             uout(i,j,k,1) = (delvx + delvy + delvz) / hbaroverm !Such that later on we can calculate uout=Lap(phase

          enddo
        enddo
      enddo

      return

      end subroutine fort_divvel

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_initcosmoax( &
                             init, i_l1,i_l2,i_l3,i_h1,i_h2,i_h3, &
                             state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                             divvel,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                             phase,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
                             didx, ivar, dx &!, maxchange &
                                 )

      use meth_params_module, only : naxvar, uaxdens, uaxre, uaxim
      use comoving_module
      use fundamental_constants_module
      use fdm_params_module

      implicit none

      integer          i_l1,i_l2,i_l3,i_h1,i_h2,i_h3
      integer          s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer          d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer          p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
      integer          didx, ivar

      double precision  init(i_l1:i_h1,i_l2:i_h2,i_l3:i_h3,ivar)
      double precision state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NAXVAR)
      double precision divvel(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,1)
      double precision phase(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,1)
      integer          i, j, k 
      

      !double precision mean_ax_dens
      double precision delta, dx(3)!, hbaroverm!, maxchange, maxchange_temp
      ! double precision weights(4) 
      ! double precision, allocatable   ::  temp(:,:,:,:)
      ! integer          s

      ! meandens = 3 * (0.7 * 100)**2 * comoving_OmAx &
      !                / ( 8 * PI * Gconst)
      !meandens = 3 * (comoving_h * 100)**2 * comoving_OmAx &
      !               / ( 8 * PI * Gconst)
      !meandens = 2.775d11 * 0.7**2 * comoving_OmAx !background density 

      !hbar/m expressed in Nyx units [Mpc km/s]
      ! hbaroverm = 0.01917152d0 / m_tt

      ! weights(:) = (/ 1.0/8.0, 1.0/16.0, 1.0/32.0, 1.0/64.0 /)
      ! s = 1

      ! allocate(temp(  i_l1:i_h1,    i_l2:i_h2,    i_l3:i_h3, 1))

      ! do n = 1,5
      !    do k = i_l3, i_h3
      !       do j = i_l2, i_h2
      !          do i = i_l1, i_h1
      !             temp(i,j,k,1) = init(i,j,k,didx) 
      !          end do
      !       end do
      !    end do
         
      !    do k = i_l3+1, i_h3-1
      !       do j = i_l2+1, i_h2-1
      !          do i = i_l1+1, i_h1-1
      !             init(i,j,k,1) = weights(1) * temp(i,j,k,1) + &         !center
      !                             weights(2) * (temp(i+s,j,k,1) + &     !sides
      !                             temp(i,j+s,k,1) + temp(i,j,k+s,1) + &
      !                             temp(i-s,j,k,1) + temp(i,j-s,k,1) + &
      !                             temp(i,j,k-s,1)) + &
      !                             weights(3) * (temp(i+s,j+s,k,1) + &   !edges
      !                             temp(i,j+s,k+s,1) + temp(i+s,j,k+s,1) + &
      !                             temp(i-s,j+s,k,1) + temp(i,j-s,k+s,1) + &
      !                             temp(i-s,j,k+s,1) + temp(i+s,j-s,k,1) + &
      !                             temp(i,j+s,k-s,1) + temp(i+s,j,k-s,1) + &
      !                             temp(i-s,j-s,k,1) + temp(i,j-s,k-s,1) + &
      !                             temp(i-s,j,k-s,1)) + &                !corners
      !                             weights(4) * (temp(i+s,j+s,k+s,1) + &
      !                             temp(i+s,j+s,k-s,1) + temp(i+s,j-s,k+s,1) + &
      !                             temp(i-s,j+s,k+s,1) + temp(i-s,j-s,k+s,1) + &
      !                             temp(i-s,j+s,k-s,1) + temp(i+s,j-s,k-s,1) + &
      !                             temp(i-s,j-s,k-s,1))
      !          end do
      !       end do
      !    end do
      ! end do
      ! deallocate(temp)

      !init has a row of ghost cells that we don't want to integrate over
      do k = i_l3+1, i_h3-1
         do j = i_l2+1, i_h2-1
            do i = i_l1+1, i_h1-1

               delta = init(i,j,k,didx)
               state(i,j,k,uaxre)   = dsqrt(delta + 1.)*dcos(phase(i,j,k,1))
               state(i,j,k,uaxim)   = dsqrt(delta + 1.)*dsin(phase(i,j,k,1))
               !state(i,j,k,UAXDENS) = meandens * (delta + 1.)
               state(i,j,k,uaxdens) = (delta + 1.)

               ! state(i,j,k,UAXRE)   = dsqrt(delta + 1.)
               ! state(i,j,k,UAXIM)   = 0.0d0
               ! state(i,j,k,UAXDENS) = meandens * (delta + 1.)

               ! state(i,j,k,UAXRE)   = divvel(i,j,k,1)
               ! state(i,j,k,UAXIM)   = phase(i,j,k,1)
               ! state(i,j,k,UAXDENS) = meandens * (state(i,j,k,UAXRE)**2+state(i,j,k,UAXIM)**2)

            enddo
         enddo
      enddo

      ! do k = i_l3+2, i_h3-2
      !    do j = i_l2+2, i_h2-2
      !       do i = i_l1+2, i_h1-2

      !          state(i,j,k,UAXIM)   = (phase(i+1,j,k,1)+phase(i-1,j,k,1)+phase(i,j+1,k,1)+phase(i,j-1,k,1)+phase(i,j,k+1,1)+phase(i,j,k-1,1)-6.d0*phase(i,j,k,1))/dx(1)/dx(1)!phase(i,j,k,1)

      !       enddo
      !    enddo
      ! enddo

      return

      end subroutine fort_initcosmoax

      
! :::
! ::: ----------------------------------------------------------------

      subroutine fort_advance_fdm_fd(time,lo,hi,&
           uin,  uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
           uout, uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
           ! grav, g_l1,g_l2,g_l3,g_h1,g_h2,g_h3, &
           phi,  p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
           delta,prob_lo,prob_hi,dt, &
           courno,a_old,a_half,a_new,verbose)

      use meth_params_module, only : naxvar, uaxdens, uaxre, uaxim
      use fdm_params_module, only : hbaroverm, ii
      use fundamental_constants_module
      use probdata_module

      implicit none

      !interface stuff
      integer          lo(3),hi(3),verbose
      integer          uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer          uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      ! integer          g_l1,g_l2,g_l3,g_h1,g_h2,g_h3
      integer          p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
      double precision  uin(  uin_l1:uin_h1,   uin_l2:uin_h2,   uin_l3:uin_h3,  NAXVAR)
      double precision uout( uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, NAXVAR)
      ! double precision grav(   g_l1:g_h1,     g_l2:g_h2,     g_l3:g_h3,   3)
      double precision  phi(   p_l1:p_h1,     p_l2:p_h2,     p_l3:p_h3)
      double precision delta(3),prob_lo(3),prob_hi(3),dt,time,courno
      double precision a_old, a_half, a_new

      !additional variables
      integer          i,j,k
      double precision invdeltasq_old(3), invdeltasq_half(3), invdeltasq_new(3)
      double precision xn,xp,xc,del,Vo,r

      !(complex) fdm field
      double complex, allocatable :: psi(:,:,:),k1(:,:,:),k2(:,:,:),k3(:,:,:),k4(:,:,:),V(:,:,:)

      allocate( psi(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3) )
      allocate(  k1(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3) )
      allocate(  k2(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3) )
      allocate(  k3(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3) )
      allocate(  k4(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3) )
      allocate(   v(  p_l1:p_h1  ,  p_l2:p_h2  ,  p_l3:p_h3  ) )

      psi = dcmplx(uin(:,:,:,uaxre),uin(:,:,:,uaxim))
      k1  = 0.d0
      k2  = 0.d0
      k3  = 0.d0
      k4  = 0.d0
      v   = 0.d0

      ! xn  = prob_hi(1)            ! Need that the center of the physical problem is at (0,0,0)!!
      ! xp  = 7.d0/8.d0*prob_hi(1)  ! Inside this radius the 'sponge' is zero.
      ! !TODO this effectively disables the sponge
      ! xp  = 1.0d10*prob_hi(1)  ! Inside this radius the 'sponge' is zero.
      ! !TODO 
      ! xc  = (xn+xp)/2.d0
      ! del = xn-xp
      ! !Vo  = 0.657d0              ! Corresponding to Vo=1 in arXiv:gr-qc/0404014v2 eq.29 when scaled to halo with rho_max=rho_cr
      ! Vo  = 0.d0

!       do k = p_l3, p_h3
!          do j = p_l2, p_h2
!             do i = p_l1, p_h1
!                r = dsqrt((i*delta(1)+prob_lo(1))**2.d0 + (j*delta(2)+prob_lo(2))**2.d0 + (k*delta(3)+prob_lo(3))**2.d0) !Distance from the physical problem center
!                if (r .lt. xp) then
!                   V(i,j,k) = 0.d0
!                else if ((r .ge. xp) .and. (r .lt. xc)) then
! !                  V(i,j,k) = ii*Vo/2.d0 * (2.d0 + dtanh((r-xc)/del) - dtanh(xc/del)) ! Cf. arXiv:gr-qc/0404014v2 eq.29 : we changed the sign so it matches the sign of phi.
!                   V(i,j,k) = ii*Vo/2.d0 * (     ((r-xp)/del)**2.d0*2.d0)
!                else
!                   V(i,j,k) = ii*Vo/2.d0 * (1.d0-((r-xn)/del)**2.d0*2.d0)
!                endif
!             enddo
!          enddo
!       enddo

      !We need this, since single scalars are not saved in checkpoint files and therefore meandens=0 after restart! 
      !if (meandens .eq. 0) then
      !   meandens = 2.775d11 * 0.7d0**2 * comoving_OmAx !3 * (100)**2 *0.7d0**2 *  comoving_OmAx / ( 8 * PI * Gconst)
      !endif

      !hbar/m expressed in Nyx units [Mpc km/s]
      !hbaroverm = 0.01917152d0 / m_tt

      do i = 1, 3
       invdeltasq_old(i)  = 1.d0 / ( a_old  * delta(i) )**2 ! differentiate w.r.t. proper distance
       invdeltasq_half(i) = 1.d0 / ( a_half * delta(i) )**2 ! differentiate w.r.t. proper distance
       invdeltasq_new(i)  = 1.d0 / ( a_new  * delta(i) )**2 ! differentiate w.r.t. proper distance
      enddo

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = uin_l3+1, uin_h3-1
         do j = uin_l2+1, uin_h2-1
            do i = uin_l1+1, uin_h1-1
               
               k1(i,j,k) = ii*hbaroverm*(6.d0*psi(i+1,j,k)&
                                        +6.d0*psi(i-1,j,k)&
                                        +6.d0*psi(i,j+1,k)&
                                        +6.d0*psi(i,j-1,k)&
                                        +6.d0*psi(i,j,k+1)&
                                        +6.d0*psi(i,j,k-1)&
                                        +3.d0*psi(i+1,j+1,k)&
                                        +3.d0*psi(i+1,j-1,k)&
                                        +3.d0*psi(i-1,j+1,k)&
                                        +3.d0*psi(i-1,j-1,k)&
                                        +3.d0*psi(i+1,j,k+1)&
                                        +3.d0*psi(i+1,j,k-1)&
                                        +3.d0*psi(i-1,j,k+1)&
                                        +3.d0*psi(i-1,j,k-1)&
                                        +3.d0*psi(i,j+1,k+1)&
                                        +3.d0*psi(i,j+1,k-1)&
                                        +3.d0*psi(i,j-1,k+1)&
                                        +3.d0*psi(i,j-1,k-1)&
                                        +2.d0*psi(i+1,j+1,k+1)&
                                        +2.d0*psi(i+1,j+1,k-1)&
                                        +2.d0*psi(i+1,j-1,k+1)&
                                        +2.d0*psi(i+1,j-1,k-1)&
                                        +2.d0*psi(i-1,j+1,k+1)&
                                        +2.d0*psi(i-1,j+1,k-1)&
                                        +2.d0*psi(i-1,j-1,k+1)&
                                        +2.d0*psi(i-1,j-1,k-1)&
                                        -88.d0*(psi(i,j,k)))&
                   *invdeltasq_old(1)/52.d0 + ii*(phi(i,j,k)+v(i,j,k))/hbaroverm*(psi(i,j,k))

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = uin_l3+1, uin_h3-1
         do j = uin_l2+1, uin_h2-1
            do i = uin_l1+1, uin_h1-1
               
               k2(i,j,k) = ii*hbaroverm*(6.d0*(psi(i+1,j,k)+k1(i+1,j,k)*dt/2.d0)&
                                        +6.d0*(psi(i-1,j,k)+k1(i-1,j,k)*dt/2.d0)&
                                        +6.d0*(psi(i,j+1,k)+k1(i,j+1,k)*dt/2.d0)&
                                        +6.d0*(psi(i,j-1,k)+k1(i,j-1,k)*dt/2.d0)&
                                        +6.d0*(psi(i,j,k+1)+k1(i,j,k+1)*dt/2.d0)&
                                        +6.d0*(psi(i,j,k-1)+k1(i,j,k-1)*dt/2.d0)&
                                        +3.d0*(psi(i+1,j+1,k)+k1(i+1,j+1,k)*dt/2.d0)&
                                        +3.d0*(psi(i+1,j-1,k)+k1(i+1,j-1,k)*dt/2.d0)&
                                        +3.d0*(psi(i-1,j+1,k)+k1(i-1,j+1,k)*dt/2.d0)&
                                        +3.d0*(psi(i-1,j-1,k)+k1(i-1,j-1,k)*dt/2.d0)&
                                        +3.d0*(psi(i+1,j,k+1)+k1(i+1,j,k+1)*dt/2.d0)&
                                        +3.d0*(psi(i+1,j,k-1)+k1(i+1,j,k-1)*dt/2.d0)&
                                        +3.d0*(psi(i-1,j,k+1)+k1(i-1,j,k+1)*dt/2.d0)&
                                        +3.d0*(psi(i-1,j,k-1)+k1(i-1,j,k-1)*dt/2.d0)&
                                        +3.d0*(psi(i,j+1,k+1)+k1(i,j+1,k+1)*dt/2.d0)&
                                        +3.d0*(psi(i,j+1,k-1)+k1(i,j+1,k-1)*dt/2.d0)&
                                        +3.d0*(psi(i,j-1,k+1)+k1(i,j-1,k+1)*dt/2.d0)&
                                        +3.d0*(psi(i,j-1,k-1)+k1(i,j-1,k-1)*dt/2.d0)&
                                        +2.d0*(psi(i+1,j+1,k+1)+k1(i+1,j+1,k+1)*dt/2.d0)&
                                        +2.d0*(psi(i+1,j+1,k-1)+k1(i+1,j+1,k-1)*dt/2.d0)&
                                        +2.d0*(psi(i+1,j-1,k+1)+k1(i+1,j-1,k+1)*dt/2.d0)&
                                        +2.d0*(psi(i+1,j-1,k-1)+k1(i+1,j-1,k-1)*dt/2.d0)&
                                        +2.d0*(psi(i-1,j+1,k+1)+k1(i-1,j+1,k+1)*dt/2.d0)&
                                        +2.d0*(psi(i-1,j+1,k-1)+k1(i-1,j+1,k-1)*dt/2.d0)&
                                        +2.d0*(psi(i-1,j-1,k+1)+k1(i-1,j-1,k+1)*dt/2.d0)&
                                        +2.d0*(psi(i-1,j-1,k-1)+k1(i-1,j-1,k-1)*dt/2.d0)&
                                        -88.d0*(psi(i,j,k)+k1(i,j,k)*dt/2.d0))&
                   *invdeltasq_half(1)/52.d0 + ii*(phi(i,j,k)+v(i,j,k))/hbaroverm*(psi(i,j,k)+k1(i,j,k)*dt/2.d0)
               
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = uin_l3+1, uin_h3-1
         do j = uin_l2+1, uin_h2-1
            do i = uin_l1+1, uin_h1-1
               
               k3(i,j,k) = ii*hbaroverm*(6.d0*(psi(i+1,j,k)+k2(i+1,j,k)*dt/2.d0)&
                                        +6.d0*(psi(i-1,j,k)+k2(i-1,j,k)*dt/2.d0)&
                                        +6.d0*(psi(i,j+1,k)+k2(i,j+1,k)*dt/2.d0)&
                                        +6.d0*(psi(i,j-1,k)+k2(i,j-1,k)*dt/2.d0)&
                                        +6.d0*(psi(i,j,k+1)+k2(i,j,k+1)*dt/2.d0)&
                                        +6.d0*(psi(i,j,k-1)+k2(i,j,k-1)*dt/2.d0)&
                                        +3.d0*(psi(i+1,j+1,k)+k2(i+1,j+1,k)*dt/2.d0)&
                                        +3.d0*(psi(i+1,j-1,k)+k2(i+1,j-1,k)*dt/2.d0)&
                                        +3.d0*(psi(i-1,j+1,k)+k2(i-1,j+1,k)*dt/2.d0)&
                                        +3.d0*(psi(i-1,j-1,k)+k2(i-1,j-1,k)*dt/2.d0)&
                                        +3.d0*(psi(i+1,j,k+1)+k2(i+1,j,k+1)*dt/2.d0)&
                                        +3.d0*(psi(i+1,j,k-1)+k2(i+1,j,k-1)*dt/2.d0)&
                                        +3.d0*(psi(i-1,j,k+1)+k2(i-1,j,k+1)*dt/2.d0)&
                                        +3.d0*(psi(i-1,j,k-1)+k2(i-1,j,k-1)*dt/2.d0)&
                                        +3.d0*(psi(i,j+1,k+1)+k2(i,j+1,k+1)*dt/2.d0)&
                                        +3.d0*(psi(i,j+1,k-1)+k2(i,j+1,k-1)*dt/2.d0)&
                                        +3.d0*(psi(i,j-1,k+1)+k2(i,j-1,k+1)*dt/2.d0)&
                                        +3.d0*(psi(i,j-1,k-1)+k2(i,j-1,k-1)*dt/2.d0)&
                                        +2.d0*(psi(i+1,j+1,k+1)+k2(i+1,j+1,k+1)*dt/2.d0)&
                                        +2.d0*(psi(i+1,j+1,k-1)+k2(i+1,j+1,k-1)*dt/2.d0)&
                                        +2.d0*(psi(i+1,j-1,k+1)+k2(i+1,j-1,k+1)*dt/2.d0)&
                                        +2.d0*(psi(i+1,j-1,k-1)+k2(i+1,j-1,k-1)*dt/2.d0)&
                                        +2.d0*(psi(i-1,j+1,k+1)+k2(i-1,j+1,k+1)*dt/2.d0)&
                                        +2.d0*(psi(i-1,j+1,k-1)+k2(i-1,j+1,k-1)*dt/2.d0)&
                                        +2.d0*(psi(i-1,j-1,k+1)+k2(i-1,j-1,k+1)*dt/2.d0)&
                                        +2.d0*(psi(i-1,j-1,k-1)+k2(i-1,j-1,k-1)*dt/2.d0)&
                                        -88.d0*(psi(i,j,k)+k2(i,j,k)*dt/2.d0))&
                   *invdeltasq_half(1)/52.d0 + ii*(phi(i,j,k)+v(i,j,k))/hbaroverm*(psi(i,j,k)+k2(i,j,k)*dt/2.d0)
               
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = uin_l3+1, uin_h3-1
         do j = uin_l2+1, uin_h2-1
            do i = uin_l1+1, uin_h1-1
               
               k4(i,j,k) = ii*hbaroverm*(6.d0*(psi(i+1,j,k)+k3(i+1,j,k)*dt)&
                                        +6.d0*(psi(i-1,j,k)+k3(i-1,j,k)*dt)&
                                        +6.d0*(psi(i,j+1,k)+k3(i,j+1,k)*dt)&
                                        +6.d0*(psi(i,j-1,k)+k3(i,j-1,k)*dt)&
                                        +6.d0*(psi(i,j,k+1)+k3(i,j,k+1)*dt)&
                                        +6.d0*(psi(i,j,k-1)+k3(i,j,k-1)*dt)&
                                        +3.d0*(psi(i+1,j+1,k)+k3(i+1,j+1,k)*dt)&
                                        +3.d0*(psi(i+1,j-1,k)+k3(i+1,j-1,k)*dt)&
                                        +3.d0*(psi(i-1,j+1,k)+k3(i-1,j+1,k)*dt)&
                                        +3.d0*(psi(i-1,j-1,k)+k3(i-1,j-1,k)*dt)&
                                        +3.d0*(psi(i+1,j,k+1)+k3(i+1,j,k+1)*dt)&
                                        +3.d0*(psi(i+1,j,k-1)+k3(i+1,j,k-1)*dt)&
                                        +3.d0*(psi(i-1,j,k+1)+k3(i-1,j,k+1)*dt)&
                                        +3.d0*(psi(i-1,j,k-1)+k3(i-1,j,k-1)*dt)&
                                        +3.d0*(psi(i,j+1,k+1)+k3(i,j+1,k+1)*dt)&
                                        +3.d0*(psi(i,j+1,k-1)+k3(i,j+1,k-1)*dt)&
                                        +3.d0*(psi(i,j-1,k+1)+k3(i,j-1,k+1)*dt)&
                                        +3.d0*(psi(i,j-1,k-1)+k3(i,j-1,k-1)*dt)&
                                        +2.d0*(psi(i+1,j+1,k+1)+k3(i+1,j+1,k+1)*dt)&
                                        +2.d0*(psi(i+1,j+1,k-1)+k3(i+1,j+1,k-1)*dt)&
                                        +2.d0*(psi(i+1,j-1,k+1)+k3(i+1,j-1,k+1)*dt)&
                                        +2.d0*(psi(i+1,j-1,k-1)+k3(i+1,j-1,k-1)*dt)&
                                        +2.d0*(psi(i-1,j+1,k+1)+k3(i-1,j+1,k+1)*dt)&
                                        +2.d0*(psi(i-1,j+1,k-1)+k3(i-1,j+1,k-1)*dt)&
                                        +2.d0*(psi(i-1,j-1,k+1)+k3(i-1,j-1,k+1)*dt)&
                                        +2.d0*(psi(i-1,j-1,k-1)+k3(i-1,j-1,k-1)*dt)&
                                        -88.d0*(psi(i,j,k)+k3(i,j,k)*dt))&
                   *invdeltasq_new(1)/52.d0 + ii*(phi(i,j,k)+v(i,j,k))/hbaroverm*(psi(i,j,k)+k3(i,j,k)*dt)
                   
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = uin_l3+1, uin_h3-1
         do j = uin_l2+1, uin_h2-1
            do i = uin_l1+1, uin_h1-1
               
               psi(i,j,k) = psi(i,j,k) + dt*(k1(i,j,k)+2.d0*k2(i,j,k)+2.d0*k3(i,j,k)+k4(i,j,k))/6.d0
               
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
         
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = uout_l3, uout_h3
        do j = uout_l2, uout_h2
          do i = uout_l1, uout_h1

             uout(i,j,k, uaxre)   = dreal(  psi(i,j,k) )  
             uout(i,j,k, uaxim)   = dimag( psi(i,j,k) )  
             !uout(i,j,k, UAXDENS) = meandens * cdabs( psi(i,j,k) )**2 !density in Nyx units
             uout(i,j,k, uaxdens) = cdabs( psi(i,j,k) )**2 !density in Nyx units

          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      deallocate(psi)
      deallocate( k1)
      deallocate( k2)
      deallocate( k3)
      deallocate( k4)
      deallocate(  v)

      end subroutine fort_advance_fdm_fd

````

