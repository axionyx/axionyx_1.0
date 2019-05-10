
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C)

      use probdata_module
      use comoving_module
      implicit none

      integer init, namlen
      integer name(namlen)
      double precision problo(3), probhi(3)

      integer untin,i

      namelist /fortin/ comoving_OmM, comoving_OmL, comoving_OmB, comoving_OmAx, comoving_h, max_num_part

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!
      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

!     Compute coordinates of grid center
      do i = 1, 3
         center(i) = (probhi(i) + problo(i)) / 2 !+ problo(i)
      end do
      end

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
!TODO_JENS: you might want to change the initialization code below in
!order to initialize the fields to sth simpler; right now this is really
!just the axion star.
      subroutine fort_initdata(level,time,lo,hi, &
                             ns, state   ,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                             na, axion,   a_l1,a_l2,a_l3,a_h1,a_h2,a_h3, &
                             nd, diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                             delta,xlo,xhi,boxsize) bind(C)
      use probdata_module
      use atomic_rates_module, only : XHYDROGEN
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT,&
                                     UFS, small_dens, TEMP_COMP, &
                                     NE_COMP, UAXDENS, UAXRE, UAXIM
      use amrex_constants_module, only : M_PI
      use axion_params_module
      use comoving_module, only : comoving_h, comoving_OmAx
      use interpolate_module

      implicit none

      integer level, ns, nd, na
      integer lo(3), hi(3)
      integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer a_l1,a_l2,a_l3,a_h1,a_h2,a_h3
      double precision xlo(3), xhi(3), time, delta(3), boxsize(3)
      double precision    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      double precision diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,nd)
      double precision    axion(a_l1:a_h1,a_l2:a_h2,a_l3:a_h3,na)

      integer i,j,k,h
      integer un,length
      double precision hubl
      double precision r,rc
      double precision rx,ry,rz
      double precision d
      double precision pi, tpi
      double precision, allocatable :: m(:), pos(:,:)
      double precision velFac, sigmaR
      double precision hbaroverm


      un = 20
      open(un,file='initial.txt',status="old",action="read")
      read(un,*) length
      read(un,*)  !Just to skip the second
      read(un,*)  !and third line

      allocate(m(1:length))
      allocate(pos(1:length,1:3))

      do i=1,length
         read(un,*) m(i), pos(i,1), pos(i,2), pos(i,3)
      end do
      close(un)
      !print *,"star is at ", pos

      ! do i=1,3
      !    pos(1,i)=center(i)
      ! enddo

      pi = 4.d0 * atan(1.d0)
      tpi = 2.0d0 * pi

      hubl = 0.7d0
      meandens = 2.775d11 * hubl**2* comoving_OmAx !background density

      rc = 1.3d0 * 0.012513007848917703d0 / (dsqrt(m_tt * hubl) * comoving_OmAx**(0.25d0))

      ! hbaroverm = 0.01917152d0 / m_tt
      velFac = 10.0d0
      sigmaR = 0.1d0

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               state(i,j,k,URHO)    = 1.5d0 * small_dens
               state(i,j,k,URHO)    = 0.0d0

               if (UMX .gt. -1) then
                  state(i,j,k,UMX:UMZ) = 0.0d0

                  ! These will both be set later in the call to init_e.
                  state(i,j,k,UEINT) = 0.d0
                  state(i,j,k,UEDEN) = 0.d0
                  diag_eos(i,j,k,TEMP_COMP) = 1000.d0
                  diag_eos(i,j,k,  NE_COMP) =    0.d0
               end if

               if (UFS .gt. -1) then
                  state(i,j,k,UFS  ) = XHYDROGEN
                  state(i,j,k,UFS+1) = (1.d0 - XHYDROGEN)
               end if
               axion(i,j,k,UAXDENS) = 0.0d0
               axion(i,j,k,UAXRE) = 0.0d0
               axion(i,j,k,UAXIM) = 0.0d0

               do h = 1,length

                   ! !!!! Previous soliton ICs
                    ! r = dsqrt((xlo(1)+(i-lo(1))*delta(1) + 0.5d0*delta(1) - pos(h,1))**2 + &
                    !           (xlo(2)+(j-lo(2))*delta(2) + 0.5d0*delta(2) - pos(h,2))**2 + &
                    !           (xlo(3)+(k-lo(3))*delta(3) + 0.5d0*delta(3) - pos(h,3))**2)
                    ! axion(i,j,k,UAXDENS) = axion(i,j,k,UAXDENS) + meandens/((1.d0+9.1d-2*(r/rc)**2.d0)**8.0d0)/length
                    ! axion(i,j,k,UAXRE)   = dsqrt(axion(i,j,k,UAXDENS))
                    ! axion(i,j,k,UAXIM)   = 0.0d0
                    !

                    rx = xlo(1)+(i-lo(1))*delta(1) + 0.5d0*delta(1) - pos(h,1)
                    ry = xlo(2)+(j-lo(2))*delta(2) + 0.5d0*delta(2) - pos(h,2)
                    rz = xlo(3)+(k-lo(3))*delta(3) + 0.5d0*delta(3) - pos(h,3)

                    ! print *, xlo(2), lo(2), pos(h,2), xhi(2), hi(2)
                    if (abs(rx) .gt. (boxsize(1)/2.0d0)) then
                       rx = boxsize(1)-abs(rx)
                    end if
                    if (abs(ry) .gt. (boxsize(2)/2.0d0)) then
                       ry = boxsize(2)-abs(ry)
                    end if
                    if (abs(rz) .gt. (boxsize(3)/2.0d0)) then
                       rz = boxsize(3)-abs(rz)
                    end if

                    r = dsqrt(rx**2 + ry**2 + rz**2)

                    ! !!!!! SINUSOIDAL DENSITY PROFILE in x direction
                    ! axion(i,j,k,UAXDENS) = 0.001*(sin(tpi*rx/xhi(1))+1.0)

                    ! !!!! GASSIAN DENSITY PROFILE
                    ! axion(i,j,k,UAXDENS) = axion(i,j,k,UAXDENS) + exp(-0.5d0*(r**2)/sigmaR/sigmaR)/sigmaR/length
                    ! !!!! GASSIAN DENSITY PROFILE - x direction
                    ! axion(i,j,k,UAXDENS) = axion(i,j,k,UAXDENS) + exp(-0.5d0*(rx**2)/sigmaR/sigmaR)/sigmaR/length

                    !!!!!! CONSTANT DENSITY profile
                    ! rx = xlo(1)+(i-lo(1))*delta(1) + 0.5d0*delta(1)
                    ! axion(i,j,k,UAXRE)   =  sin(tpi*rx/xhi(1))
                    ! axion(i,j,k,UAXIM)   =  cos(tpi*rx/xhi(1))
                    axion(i,j,k,UAXDENS) = meandens

                    !!!! SOME SOLITON-FITTING PROFILE
                    ! axion(i,j,k,UAXDENS) = axion(i,j,k,UAXDENS) + meandens/((1.d0+9.1d-2*(r/rc)**2.d0)**8.0d0)/length

                    ! !!!! CONSTANT VELOCITY KICK IN A DIAGONAL DIRECTION
                    ! axion(i,j,k,UAXRE)   = dsqrt(axion(i,j,k,UAXDENS)) * cos(velFac*(xlo(1)+(i-lo(1))*delta(1) + 0.5d0*delta(1) - pos(h,1)+xlo(2)+(j-lo(2))*delta(2) + 0.5d0*delta(2) - pos(h,2)))
                    ! axion(i,j,k,UAXIM)   = dsqrt(axion(i,j,k,UAXDENS)) * sin(velFac*(xlo(1)+(i-lo(1))*delta(1) + 0.5d0*delta(1) - pos(h,1)+xlo(2)+(j-lo(2))*delta(2) + 0.5d0*delta(2) - pos(h,2)))

                    !!!!! NO INITIAL VELOCITY
                    axion(i,j,k,UAXRE)   = dsqrt(axion(i,j,k,UAXDENS))
                    axion(i,j,k,UAXIM)   = 0.0d0

                enddo
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      deallocate(m)
      deallocate(pos)

      end subroutine fort_initdata
