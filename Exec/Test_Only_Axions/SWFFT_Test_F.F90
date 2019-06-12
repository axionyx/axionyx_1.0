#include "AMReX_LO_BCTYPES.H"

module abl_module

  use amrex_fort_module
  use amrex_error_module
  implicit none

contains

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
