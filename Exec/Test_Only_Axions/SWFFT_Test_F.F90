#include "AMReX_LO_BCTYPES.H"

module abl_module

  use amrex_fort_module
  use amrex_error_module
  implicit none

contains

  subroutine fort_init_rhs (lo, hi, rhs, rlo, rhi, domlo, domhi, problo, probhi, dx) &
       bind(c,name="fort_init_rhs")

    integer, intent(in) :: lo(3), hi(3), rlo(3), rhi(3), domlo(3), domhi(3)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in) :: problo(3), probhi(3)
    real(amrex_real), intent(in) ::  dx(3)

    integer :: i,j,k
    real(amrex_real) :: x,y,z
    real(amrex_real) :: Lx,Ly,Lz
    real(amrex_real) :: pi, fpi, tpi, fac

    Lx = probhi(1)-problo(1)
    Ly = probhi(2)-problo(2)
    Lz = probhi(3)-problo(3)

    pi = 4.d0 * atan(1.d0)
    tpi = 2.0d0 * pi
    fpi = 4.0d0 * pi
    fac = 3.0d0 * tpi**2 / (Lx**2 * Ly**2 * Lz**2)

    ! print *, "**********FORTRAN***********"
    ! print *, dx(1)
    ! print *, Lx
    do k = lo(3), hi(3)
       z = (dble(k)+0.5d0)*dx(3)/Lz

       do j = lo(2), hi(2)
          y = (dble(j)+0.5d0)*dx(2)/Ly

          do i = lo(1), hi(1)
             x = (dble(i)+0.5d0)*dx(1)/Lx

             ! rhs(i,j,k) = -fac * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                    ! &        -fac * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))

             rhs(i,j,k) = sin(tpi*x)!+sin(fpi*x)+sin(tpi*y)*sin(tpi*y)
          end do
       end do
    end do

  end subroutine fort_init_rhs

  subroutine fort_comp_asol (lo, hi, soln, slo, shi, domlo, domhi, problo, probhi, dx) &
    bind(c,name='fort_comp_asol')

    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), domlo(3), domhi(3)
    real(amrex_real), intent(inout) :: soln(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real), intent(in) :: problo(3), probhi(3)
    real(amrex_real), intent(in) :: dx(3)

    integer :: i,j,k
    real(amrex_real) :: x,y,z
    real(amrex_real) :: Lx,Ly,Lz
    real(amrex_real) :: pi, fpi, tpi

    pi = 4.d0 * atan(1.d0)
    tpi = 2.0d0 * pi
    fpi = 4.0d0 * pi

    Lx = probhi(1)-problo(1)
    Ly = probhi(1)-problo(2)
    Lz = probhi(3)-problo(3)

    do k = lo(3), hi(3)
       z = (dble(k)+0.5d0)*dx(3)
       do j = lo(2), hi(2)
          y = (dble(j)+0.5d0)*dx(2)
          do i = lo(1), hi(1)
             x = (dble(i)+0.5d0)*dx(1)

             soln(i,j,k) = 1.d0 * (sin(tpi*x/Lx) * sin(tpi*y/Ly) * sin(tpi*z/Lz))  &
                 &      + .25d0 * (sin(fpi*x/Lx) * sin(fpi*y/Ly) * sin(fpi*z/Lz))

          end do
       end do
    end do

  end subroutine fort_comp_asol

  subroutine fort_ax_fields (state, lo, hi) &
       bind(c,name="fort_ax_fields")

    use meth_params_module, only : NAXVAR, UAXDENS, UAXRE, UAXIM

    integer, intent(in) :: lo(3), hi(3)
    real(amrex_real), intent(inout) :: state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NAXVAR)

    integer :: i,j,k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             state(i,j,k,UAXDENS) = state(i,j,k,UAXRE)*state(i,j,k,UAXRE)+state(i,j,k,UAXIM)*state(i,j,k,UAXIM)
          end do
       end do
    end do

  end subroutine fort_ax_fields

  subroutine fort_kick (lo, hi, re, im, phi, hbaroverm, dt) &
       bind(c,name="fort_kick")

    integer,          intent(in)    :: lo(3), hi(3)
    real(amrex_real), intent(inout) :: re(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(amrex_real), intent(inout) :: im(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(amrex_real), intent(inout) :: phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(amrex_real), intent(in)    :: hbaroverm, dt

    integer :: i,j,k
    complex(amrex_real) :: imagi
    complex(amrex_real), allocatable :: psi(:,:,:)

    imagi = (0.0, 1.0)
    psi = dcmplx(re,im)
    psi = exp( imagi * phi / hbaroverm  * dt ) * psi
    re = real(psi)
    im = aimag(psi)

    deallocate(psi)

  end subroutine fort_kick

end module abl_module
