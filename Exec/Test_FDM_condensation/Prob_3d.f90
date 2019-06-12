
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

      dcenx = center(1)
      dceny = center(2)
      dcenz = center(3)
      dmconc = 1.0d0
      dmmass = 1.0d0
      dmscale = 1.0d0

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
                             delta,xlo,xhi,domlo,domhi) bind(C)
      use probdata_module
      use atomic_rates_module, only : XHYDROGEN
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT,&
                                     UFS, small_dens, TEMP_COMP, &
                                     NE_COMP, UAXDENS, UAXRE, UAXIM
      ! use amrex_constants_module, only : M_PI
      use fdm_params_module
      use comoving_module, only : comoving_h, comoving_OmAx
      use interpolate_module
      use fundamental_constants_module, only: Gconst, pi

      implicit none

      integer level, ns, nd, na
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer a_l1,a_l2,a_l3,a_h1,a_h2,a_h3
      double precision xlo(3), xhi(3), time, delta(3)
      double precision    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      double precision diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,nd)
      double precision    axion(a_l1:a_h1,a_l2:a_h2,a_l3:a_h3,na)

      integer i,j,k, dim, index
      integer un, domsize(3)

      do dim = 1, 3
         domsize(dim) = domhi(dim)-domlo(dim)+1
      enddo

      un = 20
      open(un,file='initial.txt',status="old",action="read")

      do index = 1,lo(3)*domsize(2)*domsize(1)
         read(un,*)
      enddo
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = lo(3), hi(3)
         do index = 1,lo(2)*domsize(1)
            read(un,*)
         enddo
         do j = lo(2), hi(2)
            do index = 1,lo(1)
               read(un,*)
            enddo
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

               read(un,*) axion(i,j,k,UAXRE), axion(i,j,k,UAXIM)
               axion(i,j,k,UAXDENS) = axion(i,j,k,UAXRE)*axion(i,j,k,UAXRE)+axion(i,j,k,UAXIM)*axion(i,j,k,UAXIM)

            enddo
            do index = 1,domhi(1)-hi(1)
               read(un,*)
            enddo
         enddo
         do index = 1,(domhi(2)-hi(2))*domsize(1)
            read(un,*)
         enddo
      enddo
      !$OMP END PARALLEL DO

      close(un)

      end subroutine fort_initdata
