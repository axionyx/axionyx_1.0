
# File Nyx\_advection\_mhd\_3d.f90

[**File List**](files.md) **>** [**MHD**](dir_a5db59ee0cc93a408ad0433ba32613c6.md) **>** [**Nyx\_advection\_mhd\_3d.f90**](Nyx__advection__mhd__3d_8f90.md)

[Go to the documentation of this file.](Nyx__advection__mhd__3d_8f90.md) 


````cpp
! :::
! ::: ----------------------------------------------------------------
! :::

!fort_make_mhd_sources generates the magnetohydrodynamical fluxes and combines them with Source 
!terms (TODO) 
!Currently Computes Fluxes combines them into "hydro_src" and computes Electric Fields 
  subroutine fort_make_mhd_sources(time, lo, hi, &
           uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
           bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
           byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
           bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
           bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
           byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
           bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
           ugdnvx,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
           ugdnvy,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
           ugdnvz,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
           src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
           hydro_src, hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3, &
           divu_cc, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
           grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
           delta , dt, &
           flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
           flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
           flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
           Ex,ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
           Ey,ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
           Ez,ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
           courno, a_old,a_new,print_fortran_warnings) &
           bind(C, name="fort_make_mhd_sources")
        
 !--------------------- Dependencies ------------------------------------------------
      use amrex_fort_module, only : rt => amrex_real
      use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
      use ct_upwind, only : corner_transport, checkisnan
      use mhd_plm_module, only : plm
      use meth_params_module!, only : QVAR, NTHERM, NHYP, normalize_species, NVAR, URHO, UEDEN
      use enforce_module, only : enforce_nonnegative_species
      use bl_constants_module

      implicit none

!-------------------- Variables -----------------------------------------------------

      integer lo(3),hi(3),print_fortran_warnings,do_grav
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
      integer byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
      integer bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
      integer bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3
      integer byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3
      integer bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3
      integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
      integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
      integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      integer ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3
      integer ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3
      integer ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
      integer d_l1, d_l2, d_l3, d_h1, d_h2, d_h3

      real(rt)  uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3,  NVAR)
      real(rt)  bxin(bxin_l1:bxin_h1, bxin_l2:bxin_h2, bxin_l3:bxin_h3)
      real(rt)  bxout(bxout_l1:bxout_h1, bxout_l2:bxout_h2, bxout_l3:bxout_h3)
      real(rt)  byin(byin_l1:byin_h1, byin_l2:byin_h2, byin_l3:byin_h3)
      real(rt)  byout(byout_l1:byout_h1, byout_l2:byout_h2, byout_l3:byout_h3)
      real(rt)  bzin(bzin_l1:bzin_h1, bzin_l2:bzin_h2, bzin_l3:bzin_h3)
      real(rt)  bzout(bzout_l1:bzout_h1, bzout_l2:bzout_h2, bzout_l3:bzout_h3)
      real(rt)  src(src_l1:src_h1, src_l2:src_h2, src_l3:src_h3, NTHERM)
      real(rt)  hydro_src(hsrc_l1:hsrc_h1, hsrc_l2:hsrc_h2, hsrc_l3:hsrc_h3, NVAR)
      real(rt)  divu_cc(d_l1:d_h1, d_l2:d_h2, d_l3:d_h3)
      real(rt)  ugdnvx(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
      real(rt)  ugdnvy(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
      real(rt)  ugdnvz(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
      real(rt)  grav( gv_l1:gv_h1, gv_l2:gv_h2, gv_l3:gv_h3, 3)

      real(rt), intent(inout) ::  flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,NVAR)
      real(rt), intent(inout) ::  flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,NVAR)
      real(rt), intent(inout) ::  flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,NVAR)

      real(rt), intent(inout) ::  Ex(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)
      real(rt), intent(inout) ::  Ey(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)
      real(rt), intent(inout) ::  Ez(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)

      real(rt)  delta(3),dt,time, courno
      real(rt)  a_old, a_new

      integer flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
      integer flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
      integer flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

      integer extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3
      integer eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3
      integer eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3

      ! Automatic arrays for workspace
      real(rt), pointer :: q(:,:,:,:)
      real(rt), pointer :: cx(:,:,:)
      real(rt), pointer :: cy(:,:,:)
      real(rt), pointer :: cz(:,:,:)
      real(rt), pointer :: csml(:,:,:)
      real(rt), pointer :: flatn(:,:,:)
      real(rt), pointer :: srcQ(:,:,:,:)

      real(rt), allocatable :: flxx(:,:,:,:)
      real(rt), allocatable :: flxy(:,:,:,:)
      real(rt), allocatable :: flxz(:,:,:,:)

      real(rt), allocatable :: Extemp(:,:,:)
      real(rt), allocatable :: Eytemp(:,:,:)
      real(rt), allocatable :: Eztemp(:,:,:)

      real(rt), allocatable :: qp(:,:,:,:,:)
      real(rt), allocatable :: qm(:,:,:,:,:)

      real(rt) dx,dy,dz
      integer ngq,ngf
      integer q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3
      integer :: i,j,k
      real(rt), pointer :: bcc(:,:,:,:)

      call amrex_allocate(bcc, lo-nhyp, hi+nhyp, 3)
   
      ngq = nhyp
      ngf = 1
      q_l1 = lo(1)-nhyp
      q_l2 = lo(2)-nhyp
      q_l3 = lo(3)-nhyp
      q_h1 = hi(1)+nhyp
      q_h2 = hi(2)+nhyp
      q_h3 = hi(3)+nhyp

      srcq_l1 = lo(1)-1
      srcq_l2 = lo(2)-1
      srcq_l3 = lo(3)-1
      srcq_h1 = hi(1)+1
      srcq_h2 = hi(2)+1
      srcq_h3 = hi(3)+1

      call amrex_allocate(     q, lo-nhyp, hi+nhyp, qvar)
      call amrex_allocate( flatn, lo-nhyp, hi+nhyp      )
      call amrex_allocate(    cx, lo-nhyp, hi+nhyp      )
      call amrex_allocate(    cy, lo-nhyp, hi+nhyp      )
      call amrex_allocate(    cz, lo-nhyp, hi+nhyp      )
      call amrex_allocate(  csml, lo-nhyp, hi+nhyp      )
      call amrex_allocate(  srcq, lo-1, hi+1, qvar)

      flxx_l1 = lo(1)-3
      flxx_l2 = lo(2)-3
      flxx_l3 = lo(3)-3
      flxx_h1 = hi(1)+4
      flxx_h2 = hi(2)+3
      flxx_h3 = hi(3)+3

      flxy_l1 = lo(1)-3
      flxy_l2 = lo(2)-3
      flxy_l3 = lo(3)-3
      flxy_h1 = hi(1)+3
      flxy_h2 = hi(2)+4
      flxy_h3 = hi(3)+3

      flxz_l1 = lo(1)-3
      flxz_l2 = lo(2)-3
      flxz_l3 = lo(3)-3
      flxz_h1 = hi(1)+3
      flxz_h2 = hi(2)+3
      flxz_h3 = hi(3)+4

      extemp_l1 = lo(1)-3
      extemp_l2 = lo(2)-3
      extemp_l3 = lo(3)-3
      extemp_h1 = hi(1)+3
      extemp_h2 = hi(2)+4
      extemp_h3 = hi(3)+4

      eytemp_l1 = lo(1)-3
      eytemp_l2 = lo(2)-3
      eytemp_l3 = lo(3)-3
      eytemp_h1 = hi(1)+4
      eytemp_h2 = hi(2)+3
      eytemp_h3 = hi(3)+4

      eztemp_l1 = lo(1)-3
      eztemp_l2 = lo(2)-3
      eztemp_l3 = lo(3)-3
      eztemp_h1 = hi(1)+4
      eztemp_h2 = hi(2)+4
      eztemp_h3 = hi(3)+3

      allocate(flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,qvar))
      allocate(flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,qvar))
      allocate(flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,qvar))

      allocate(extemp(extemp_l1:extemp_h1,extemp_l2:extemp_h2,extemp_l3:extemp_h3))
      allocate(eytemp(eytemp_l1:eytemp_h1,eytemp_l2:eytemp_h2,eytemp_l3:eytemp_h3))
      allocate(eztemp(eztemp_l1:eztemp_h1,eztemp_l2:eztemp_h2,eztemp_l3:eztemp_h3))

      allocate(  qp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,qvar, 3))
      allocate(  qm(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,qvar, 3))

      q = 0.d0

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

!Calculate Primitives based on conservatives
      call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3,&
                bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
                byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
                bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
                q , cx , cy, cz , csml, flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                src,  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                srcq, srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                grav,gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                courno, dx,dy,dz,dt,ngq,ngf,a_old,a_new)

!Generate Cell Centered Magnetic Fields for Energy Correction 
!gen_bcc creates the cell centered Magnetic Fields at the current time step, to be updated in flux_ener_corr 
      call gen_bcc(lo, hi, bcc, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                   bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
                   byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
                   bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3)

!Interpolate Cell centered values to faces
      call plm(lo, hi, q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3,&
               bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
               byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
               bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
               qp, qm, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, dx, dy, dz, dt, a_old)
!      qm(:,:,:,:,1) = q
!      qp(:,:,:,:,1) = q
!      qm(:,:,:,:,2) = q
!      qp(:,:,:,:,2) = q
!      qm(:,:,:,:,3) = q
!      qp(:,:,:,:,3) = q

      flxx = 0.d0
      flxy = 0.d0
      flxz = 0.d0

!Corner Couple and find the correct fluxes + electric fields
      call corner_transport( q, qm, qp, q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3, &
                 flxx,flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                 flxy,flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                 flxz,flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                 extemp, extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3, &
                 eytemp, eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3, &
                 eztemp, eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3, &
                 dx , dy, dz, dt)

     flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,urho:ueden) = &
     flxx(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,urho:ueden)

     flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,urho:ueden) = &
     flxy(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,urho:ueden)

     flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,urho:ueden) = &
     flxz(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,urho:ueden)

!Magnetic Update
     call magup(bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
                byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
                bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
                bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
                byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
                bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
                src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
                extemp, extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3, &
                eytemp, eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3, &
                eztemp, eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3, &
                lo, hi, dx, dy, dz, dt, a_old, a_new)

     ex(ex_l1:ex_h1,ex_l2:ex_h2, ex_l3:ex_h3) = extemp(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)
     ey(ey_l1:ey_h1,ey_l2:ey_h2, ey_l3:ey_h3) = eytemp(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)
     ez(ez_l1:ez_h1,ez_l2:ez_h2, ez_l3:ez_h3) = eztemp(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)

!Should be the end of Make MHD Sources, export the Electric Field and combine the Fluxes into hydro_sources
     call flux_combo(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                      hydro_src, hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3, &
                      flux1, flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3, &
                      flux2, flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3, &
                      flux3, flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3, &
                      divu_cc, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, lo, hi, & 
                      dx, dy, dz, dt, a_old, a_new)                  

! Removes the non-solenoidal contribution to the total energy flux
    call flux_ener_corr(hydro_src, hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3, & 
                        bcc, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, & 
                        flxx,flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                        flxy,flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                        flxz,flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                        bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
                        byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
                        bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
                        lo, hi, dx, dy, dz, dt, a_old, a_new)
         
! We are done with these here so can go ahead and free up the space
     call amrex_deallocate(q)
     call amrex_deallocate(flatn)
     call amrex_deallocate(cx)
     call amrex_deallocate(cy)
     call amrex_deallocate(cz)
     call amrex_deallocate(csml)
     call amrex_deallocate(bcc)
!    call bl_deallocate(div)
     call amrex_deallocate(srcq)
!    call bl_deallocate(pdivu)

     deallocate(qm)
     deallocate(qp)

     deallocate(flxx,flxy,flxz)
     deallocate(extemp,eytemp,eztemp)

  end subroutine fort_make_mhd_sources

!Updates Hydrodynamical Variables with corrected "hydro_sources" notably fluxes
  subroutine fort_update_mhd_state (lo, hi, uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                  uout      , uout_l1 , uout_l2 , uout_l3 , uout_h1 , uout_h2 , uout_h3 , &
                  bxout     , bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
                  byout     , byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
                  bzout     , bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
                  src       , src_l1  , src_l2  , src_l3  , src_h1  , src_h2  , src_h3  , &
                  hydro_src , hsrc_l1 , hsrc_l2 , hsrc_l3 , hsrc_h1 , hsrc_h2 , hsrc_h3 , & 
                  divu_cc   , d_l1    , d_l2    , d_l3    , d_h1    , d_h2    , d_h3    , &
                  dt        , dx      , dy      , dz      ,a_old    , a_new   , print_fortran_warnings) &
                  bind(c, name='fort_update_mhd_state')

!--------------------- Dependencies ------------------------------------------------
      use amrex_fort_module, only : rt => amrex_real
      use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
      use meth_params_module!, only : QVAR, NTHERM, NHYP, normalize_species, NVAR, URHO, UEDEN
      use enforce_module, only : enforce_nonnegative_species
      use bl_constants_module
      implicit none

!-------------------- Variables -----------------------------------------------------

      integer lo(3),hi(3),print_fortran_warnings
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3
      integer d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
      integer bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3
      integer byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3
      integer bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3


      real(rt)  uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3,  NVAR)
      real(rt)  uout(uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, NVAR)
      real(rt)  src(src_l1:src_h1, src_l2:src_h2, src_l3:src_h3, NTHERM)
      real(rt)  bxout(bxout_l1:bxout_h1, bxout_l2:bxout_h2, bxout_l3:bxout_h3)
      real(rt)  byout(byout_l1:byout_h1, byout_l2:byout_h2, byout_l3:byout_h3)
      real(rt)  bzout(bzout_l1:bzout_h1, bzout_l2:bzout_h2, bzout_l3:bzout_h3)

      real(rt)  hydro_src(hsrc_l1:hsrc_h1, hsrc_l2:hsrc_h2, hsrc_l3:hsrc_h3, NVAR)
      real(rt)  divu_cc(d_l1:d_h1, d_l2:d_h2, d_l3:d_h3)
 
      real(rt)  dt, a_old, a_new, dx, dy, dz, i, j, k 
    
!Conservative update
      call consup(uin,  uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                  uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                  bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
                  byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
                  bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
                  hydro_src, hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3, &
                  lo ,hi , dt ,a_old ,a_new)
      
              
! Enforce the density >= small_dens.  Make sure we do this immediately after updates.
    call enforce_minimum_density(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                                 uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                 bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
                                 byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
                                 bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
                                 lo,hi,print_fortran_warnings)
!      if (do_grav .gt. 0)  then
!          call add_grav_source(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
!                               uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
!                               grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
!                               lo,hi,dx,dy,dz,dt,a_old,a_new,e_added,ke_added)
!      endif
!     Enforce species >= 0
!     call enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
!                                      uout_h1,uout_h2,uout_h3,lo,hi,0)

!     Re-normalize the species
!     if (normalize_species .eq. 1) then
!         call normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
!                               lo,hi)
!     end if
  end subroutine fort_update_mhd_state       


  subroutine ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3,&
                         bx, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
                         by, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
                         bz, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
                         q,cx, cy, cz,csml,flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                         src,  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                         srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                         grav,gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                         courno, dx,dy,dz,dt,ngp,ngf,a_old,a_new)
      !
      !     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
      !     if use_flattening=1.  Declared dimensions of q,c,csml,flatn are given
      !     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
      !     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
      !     routine that computes flatn).
      !
      use amrex_fort_module, only : rt => amrex_real
      use network, only : nspec, naux
      use eos_params_module
      use eos_module
!      use flatten_module
      use bl_constants_module
      use meth_params_module

      implicit none

      real(rt), parameter:: small = 1.d-8

      integer lo(3), hi(3)
      integer  uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
      integer bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
      integer byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
      integer bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
      integer    q_l1,   q_l2,   q_l3,   q_h1,   q_h2,   q_h3
      integer   gv_l1,  gv_l2,  gv_l3,  gv_h1,  gv_h2,  gv_h3
      integer  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
      integer srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3

      real(rt) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NTHERM)
      real(rt) :: bx(bxin_l1:bxin_h1, bxin_l2:bxin_h2, bxin_l3:bxin_h3)
      real(rt) :: by(byin_l1:byin_h1, byin_l2:byin_h2, byin_l3:byin_h3)
      real(rt) :: bz(bzin_l1:bzin_h1, bzin_l2:bzin_h2, bzin_l3:bzin_h3)

      real(rt) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR) !Contains Cell Centered Mag Field
      real(rt) :: cx(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: cy(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: cz(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: csml(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: src( src_l1: src_h1, src_l2: src_h2, src_l3: src_h3,NTHERM)
      real(rt) :: srcQ(srcq_l1:srcq_h1,srcq_l2:srcq_h2,srcq_l3:srcq_h3,QVAR)
      real(rt) :: grav( gv_l1: gv_h1, gv_l2: gv_h2, gv_l3: gv_h3,3)
      real(rt) :: dx, dy, dz, dt, courno, a_old, a_new
      real(rt) :: dpdr, dpde

      integer          :: i, j, k
      integer          :: ngp, ngf, loq(3), hiq(3)
      integer          :: n, nq
      integer          :: iadv, ispec
      real(rt) :: courx, coury, courz, courmx, courmy, courmz, cad
      real(rt) :: a_half, a_dot, rhoInv
      real(rt) :: dtdxaold, dtdyaold, dtdzaold, small_pres_over_dens

      do i=1,3
         loq(i) = lo(i)-ngp
         hiq(i) = hi(i)+ngp
      enddo
      !
      ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
      ! The temperature is used as an initial guess for the eos call and will be overwritten.
      !
      !Calculate Cell Centered Magnetic Field x

      do k = loq(3), hiq(3)
         do j = loq(2), hiq(2)
            do i = loq(1), hiq(1)
              q(i,j,k,qmagx) = 0.5d0*(bx(i+1,j,k) + bx(i,j,k))
            end do
         end do
      end do

      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
             q(i,j,k,qmagy) = 0.5d0*(by(i,j+1,k) + by(i,j,k))
            end do
         end do
      end do

      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
             q(i,j,k,qmagz) = 0.5d0*(bz(i,j,k+1) + bz(i,j,k))
            end do
         end do
      end do
      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               if (uin(i,j,k,urho) .le. zero) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> Error: Nyx_advection_3d::ctoprim ',i,j,k
                  print *,'>>> ... negative density ',uin(i,j,k,urho)
                  call bl_error("Error:: Nyx_advection_3d.f90 :: ctoprim")
               end if

               rhoinv = one/uin(i,j,k,urho)

               q(i,j,k,qrho) = uin(i,j,k,urho)
               q(i,j,k,qu)   = uin(i,j,k,umx)*rhoinv
               q(i,j,k,qv)   = uin(i,j,k,umy)*rhoinv
               q(i,j,k,qw)   = uin(i,j,k,umz)*rhoinv

               ! Convert "rho e" to "e"
               q(i,j,k,qreint ) = uin(i,j,k,ueint)*rhoinv
            enddo
         enddo
      enddo

      ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
      do iadv = 1, nadv
         n = ufa + iadv - 1
         nq = qfa + iadv - 1
         do k = loq(3),hiq(3)
            do j = loq(2),hiq(2)
               do i = loq(1),hiq(1)
                  q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,qrho)
               enddo
            enddo
         enddo
      enddo

      ! Load chemical species and aux. variables, c, into q, assuming they arrived in uin as rho.c
!      if (UFS .gt. 0) then
!         do ispec = 1, nspec+naux
!            n  = UFS + ispec - 1
!            nq = QFS + ispec - 1
!            do k = loq(3),hiq(3)
!               do j = loq(2),hiq(2)
!                  do i = loq(1),hiq(1)
!                     q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
!                  enddo
!               enddo
!            enddo
!         enddo
!      end if ! UFS > 0

      small_pres_over_dens = small_pres / small_dens


      ! Get p, T, c, csml using q state
      do k = loq(3), hiq(3)
         do j = loq(2), hiq(2)
            do i = loq(1), hiq(1)

               ! If necessary, reset the energy using small_temp
               if (q(i,j,k,qreint) .lt. zero) then

!                 HACK HACK HACK 
!                 call nyx_eos_given_RT(q(i,j,k,QREINT),q(i,j,k,QPRES),q(i,j,k,QRHO), &
!                                       small_temp,diag_eos(i,j,k,NE_COMP),a_old)

                  if (q(i,j,k,qreint) .lt. zero) then
                     !
                     ! A critical region since we usually can't write from threads.
                     !
                     print *,'   '
                     print *,'>>> Error: Nyx_advection_3d::ctoprim ',i,j,k
                     print *,'>>> ... new e from eos_given_RT call is negative ',q(i,j,k,qreint)
                     print *,'    '
                     call bl_error("Error:: Nyx_advection_3d.f90 :: ctoprim")
                  end if
               end if

               ! Define the magneto-accoustic speed from the EOS
               cad = q(i,j,k,qmagx)!(q(i,j,k,QMAGX)**2)/q(i,j,k,QRHO)
               call nyx_eos_soundspeed(cx(i,j,k), q(i,j,k,qrho), q(i,j,k,qreint), &
               q(i,j,k,qmagx), q(i,j,k,qmagy), q(i,j,k,qmagz), cad)

               cad = q(i,j,k,qmagy)!(q(i,j,k,QMAGY)**2)/q(i,j,k,QRHO)
               call nyx_eos_soundspeed(cy(i,j,k), q(i,j,k,qrho), q(i,j,k,qreint), &
               q(i,j,k,qmagx), q(i,j,k,qmagy), q(i,j,k,qmagz), cad)

               cad = q(i,j,k,qmagz)!(q(i,j,k,QMAGZ)**2)/q(i,j,k,QRHO)
               call nyx_eos_soundspeed(cz(i,j,k), q(i,j,k,qrho), q(i,j,k,qreint), &
               q(i,j,k,qmagx), q(i,j,k,qmagy), q(i,j,k,qmagz), cad)

               ! Set csmal based on small_pres and small_dens
               csml(i,j,k) = sqrt(gamma_const * small_pres_over_dens)

               ! Convert "e" back to "rho e"
               q(i,j,k,qreint) = q(i,j,k,qreint)*q(i,j,k,qrho)

               ! Pressure = (gamma - 1) * rho * e
               q(i,j,k,qpres) = gamma_minus_1 * q(i,j,k,qreint)
            end do
         end do
      end do

      a_half = half * (a_old + a_new)
      a_dot   = (a_new - a_old) / dt

      ! Make sure these are initialized to zero.
      srcq = zero

      ! NOTE - WE ASSUME HERE THAT src(i,j,k,URHO) = 0. --
      !        IF NOT THEN THE FORMULAE BELOW ARE INCOMPLETE.

      ! compute srcQ terms
  !    do k = lo(3)-1, hi(3)+1
  !       do j = lo(2)-1, hi(2)+1
  !          do i = lo(1)-1, hi(1)+1

  !             rhoInv = ONE/q(i,j,k,QRHO)

  !            srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
  !             srcQ(i,j,k,QU    ) = src(i,j,k,UMX) * rhoInv - a_dot * q(i,j,k,QU) + &
  !                                  grav(i,j,k,1)
  !             srcQ(i,j,k,QV    ) = src(i,j,k,UMY) * rhoInv - a_dot * q(i,j,k,QV) + &
  !                                  grav(i,j,k,2)
  !             srcQ(i,j,k,QW    ) = src(i,j,k,UMZ) * rhoInv - a_dot * q(i,j,k,QW) + &
  !                                  grav(i,j,k,3)
  !             srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) - &
  !                                                     q(i,j,k,QV)*src(i,j,k,UMY) - &
  !                                                     q(i,j,k,QW)*src(i,j,k,UMZ) - &
  !                                                     a_dot * THREE * gamma_minus_1 * q(i,j,k,QREINT)

   !            dpde = gamma_minus_1 * q(i,j,k,QRHO)
   !            dpdr = gamma_minus_1 * q(i,j,k,QREINT)/q(i,j,k,QRHO)
   !            srcQ(i,j,k,QPRES ) = dpde * srcQ(i,j,k,QREINT) * rhoInv &
   !                               + dpdr * srcQ(i,j,k,QRHO)

   !            if (UFS .gt. 0) then
   !               do ispec = 1,nspec+naux
   !                  srcQ(i,j,k,QFS+ispec-1) = src(i,j,k,UFS+ispec-1)*rhoInv
   !               enddo
   !            end if ! UFS > 0

   !            do iadv = 1,nadv
   !               srcQ(i,j,k,QFA+iadv-1) = src(i,j,k,UFA+iadv-1)*rhoInv
   !            enddo

   !         enddo
   !      enddo
   !   enddo

      ! Compute running max of Courant number over grids
      courmx = courno
      courmy = courno
      courmz = courno

      dtdxaold = dt / dx / a_old
      dtdyaold = dt / dy / a_old
      dtdzaold = dt / dz / a_old

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               courx = ( cx(i,j,k)+abs(q(i,j,k,qu)) ) * dtdxaold
               coury = ( cy(i,j,k)+abs(q(i,j,k,qv)) ) * dtdyaold
               courz = ( cz(i,j,k)+abs(q(i,j,k,qw)) ) * dtdzaold

               courmx = max( courmx, courx )
               courmy = max( courmy, coury )
               courmz = max( courmz, courz )

            enddo
         enddo
      enddo
      courno = max( courmx, courmy, courmz )

      end subroutine ctoprim
! :::
! ::: ========================== Conservative Update ===============================================================
! ::: 

      subroutine consup(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                        uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                        bx, bx_l1, bx_l2, bx_l3, bx_h1, bx_h2, bx_h3, &
                        by, by_l1, by_l2, by_l3, by_h1, by_h2, by_h3, &
                        bz, bz_l1, bz_l2, bz_l3, bz_h1, bz_h2, bz_h3, &
                        hydro_src ,  hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3, &
                        lo,hi,dt,a_old,a_new)

       use amrex_fort_module, only : rt => amrex_real
       use meth_params_module, only : nvar, urho, umx, umy, umz, ueint, ueden

      implicit none

      integer,  intent(in)    :: uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
      integer,  intent(in)    :: bx_l1,  bx_l2,  bx_l3,  bx_h1,  bx_h2,  bx_h3
      integer,  intent(in)    :: by_l1,  by_l2,  by_l3,  by_h1,  by_h2,  by_h3 
      integer,  intent(in)    :: bz_l1,  bz_l2,  bz_l3,  bz_h1,  bz_h2,  bz_h3
      integer,  intent(in)    :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer,  intent(in)    :: hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3
      integer,  intent(in)    :: lo(3), hi(3)

      real(rt), intent(in)    :: uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3, NVAR)
      real(rt), intent(in)    :: hydro_src(hsrc_l1:hsrc_h1,hsrc_l2:hsrc_h2,hsrc_l3:hsrc_h3, NVAR)
      real(rt), intent(in)    :: bx(bx_l1:bx_h1,bx_l2:bx_h2,bx_l3:bx_h3)
      real(rt), intent(in)    :: by(by_l1:by_h1,by_l2:by_h2,by_l3:by_h3)
      real(rt), intent(in)    :: bz(bz_l1:bz_h1,bz_l2:bz_h2,bz_l3:bz_h3)
      real(rt), intent(in)    :: dt,a_old, a_new 
      real(rt), intent(out)   :: uout(uout_l1:uout_h1,uout_l2:uout_h2, uout_l3:uout_h3,NVAR)
     

      integer                 :: i, j, k, n
      real(rt)                :: u, v, w, rhoinv, bcx, bcy, bcz
      do n = 1, nvar
        do k = lo(3), hi(3)
          do j = lo(2), hi(2)
            do i = lo(1), hi(1)
              uout(i,j,k,n) = uin(i,j,k,n) + dt*(hydro_src(i,j,k,n))
            enddo
          enddo
        enddo
      enddo
      
!Generate UEINT from UEDEN, other sources TODO 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
              rhoinv = 1.d0/uout(i,j,k,urho)
              u = uout(i,j,k,umx)*rhoinv
              v = uout(i,j,k,umy)*rhoinv
              w = uout(i,j,k,umz)*rhoinv 
              bcx = 0.5d0*(bx(i,j,k) + bx(i+1,j,k))
              bcy = 0.5d0*(by(i,j,k) + by(i,j+1,k))
              bcz = 0.5d0*(bz(i,j,k) + bz(i,j,k+1)) 
              uout(i,j,k,ueint) = uout(i,j,k,ueden) - 0.5d0*(u*u + v*v + w*w)*uout(i,j,k,urho) - 0.5d0*(bcx*bcx + bcy*bcy + bcz*bcz)
            enddo
         enddo
      enddo
      end subroutine consup

! :::
! ::: ========================== Magnetic Update ===============================================================
! ::: 

  subroutine magup(bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
                   byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
                   bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
                   bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
                   byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
                   bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
                   src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
                   Ex,ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                   Ey,ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                   Ez,ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                   lo, hi, dx, dy, dz, dt, a_old, a_new)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module!, only : QVAR, NVAR, UEINT

  implicit none
  
    integer, intent(in)   :: bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
    integer, intent(in)   :: byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
    integer, intent(in)   :: bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
    integer, intent(in)   :: bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3
    integer, intent(in)   :: byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3
    integer, intent(in)   :: bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3
    integer, intent(in)   :: src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3
    integer, intent(in)   :: ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3
    integer, intent(in)   :: ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3
    integer, intent(in)   :: ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3
    integer, intent(in)   :: lo(3), hi(3)

    real(rt), intent(in)  :: bxin(bxin_l1:bxin_h1, bxin_l2:bxin_h2, bxin_l3:bxin_h3)
    real(rt), intent(in)  :: byin(byin_l1:byin_h1, byin_l2:byin_h2, byin_l3:byin_h3)
    real(rt), intent(in)  :: bzin(bzin_l1:bzin_h1, bzin_l2:bzin_h2, bzin_l3:bzin_h3)
    real(rt), intent(in)  :: src(src_l1:src_h1, src_l2:src_h2, src_l3:src_h3, QVAR)

    real(rt), intent(in)  ::  Ex(ex_l1:ex_h1,ex_l2:ex_h2, ex_l3:ex_h3)
    real(rt), intent(in)  ::  Ey(ey_l1:ey_h1,ey_l2:ey_h2, ey_l3:ey_h3)
    real(rt), intent(in)  ::  Ez(ez_l1:ez_h1,ez_l2:ez_h2, ez_l3:ez_h3)

    real(rt), intent(in)  :: dx, dy, dz, dt, a_old, a_new 

    real(rt), intent(out) :: bxout(bxout_l1:bxout_h1, bxout_l2:bxout_h2, bxout_l3:bxout_h3)
    real(rt), intent(out) :: byout(byout_l1:byout_h1, byout_l2:byout_h2, byout_l3:byout_h3)
    real(rt), intent(out) :: bzout(bzout_l1:bzout_h1, bzout_l2:bzout_h2, bzout_l3:bzout_h3)

    real(rt)              :: dtdxinv, dtdyinv ,dtdzinv 
    integer               :: i, j, k

    dtdxinv = dt/dx
    dtdyinv = dt/dy 
    dtdzinv = dt/dz
  !***** TO DO ***** SOURCES
  !-------------------------------- bx --------------------------------------------------
    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
          bxout(i,j,k) = bxin(i,j,k) + dtdzinv*(ey(i,j,k+1) - ey(i,j,k)) - dtdyinv*(ez(i,j+1,k) - ez(i,j,k))
        enddo
      enddo
    enddo

  !------------------------------- by --------------------------------------------------
  do k = lo(3), hi(3)
    do j = lo(2), hi(2)+1
      do i = lo(1), hi(1)
        byout(i,j,k) = byin(i,j,k) + dtdxinv*(ez(i+1,j,k) - ez(i,j,k)) - dtdzinv*(ex(i,j,k+1) - ex(i,j,k))
      enddo
    enddo
  enddo
  !------------------------------- bz --------------------------------------------------
  do k = lo(3), hi(3)+1
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        bzout(i,j,k) = bzin(i,j,k) + dtdyinv*(ex(i,j+1,k) - ex(i,j,k)) - dtdxinv*(ey(i+1,j,k) - ey(i,j,k))
      enddo
    enddo
  enddo
  end subroutine magup

! :::
! ::: ========================== Energy Correction ===========================================================
! ::: 

  subroutine flux_ener_corr(hydro_src, hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3, & 
                      bcc , bcc_l1, bcc_l2, bcc_l3, bcc_h1, bcc_h2, bcc_h3, &
                      flxx, flxx_l1, flxx_l2, flxx_l3, flxx_h1, flxx_h2, flxx_h3, &
                      flxy, flxy_l1, flxy_l2, flxy_l3, flxy_h1, flxy_h2, flxy_h3, &
                      flxz, flxz_l1, flxz_l2, flxz_l3, flxz_h1, flxz_h2, flxz_h3, &
                      bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
                      byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
                      bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
                      lo, hi, dx, dy, dz, dt, a_old, a_new)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module

  implicit none
  
  integer,  intent(in   )   :: hsrc_l1 , hsrc_l2 , hsrc_l3 , hsrc_h1 , hsrc_h2 , hsrc_h3 
  integer,  intent(in   )   :: flxx_l1 , flxx_l2 , flxx_l3 , flxx_h1 , flxx_h2 , flxx_h3
  integer,  intent(in   )   :: flxy_l1 , flxy_l2 , flxy_l3 , flxy_h1 , flxy_h2 , flxy_h3
  integer,  intent(in   )   :: flxz_l1 , flxz_l2 , flxz_l3 , flxz_h1 , flxz_h2 , flxz_h3
  integer,  intent(in   )   :: bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3
  integer,  intent(in   )   :: byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3
  integer,  intent(in   )   :: bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3
  integer,  intent(in   )   :: bcc_l1, bcc_l2, bcc_l3, bcc_h1, bcc_h2, bcc_h3
  integer,  intent(in   )   :: lo(3), hi(3)

  real(rt), intent(inout) :: hydro_src(hsrc_l1:hsrc_h1,hsrc_l2:hsrc_h2, hsrc_l3:hsrc_h3, NVAR)
  real(rt), intent(in   ) :: flxx(flxx_l1:flxx_h1, flxx_l2:flxx_h2, flxx_l3:flxx_h3, QVAR)
  real(rt), intent(in   ) :: flxy(flxy_l1:flxy_h1, flxy_l2:flxy_h2, flxy_l3:flxy_h3, QVAR)
  real(rt), intent(in   ) :: flxz(flxz_l1:flxz_h1, flxz_l2:flxz_h2, flxz_l3:flxz_h3, QVAR)
    real(rt), intent(in   ) :: bxout(bxout_l1:bxout_h1, bxout_l2:bxout_h2, bxout_l3:bxout_h3)
    real(rt), intent(in   ) :: byout(byout_l1:byout_h1, byout_l2:byout_h2, byout_l3:byout_h3)
    real(rt), intent(in   ) :: bzout(bzout_l1:bzout_h1, bzout_l2:bzout_h2, bzout_l3:bzout_h3)
  real(rt), intent(inout) :: bcc(bcc_l1:bcc_h1,bcc_l2:bcc_h2, bcc_l3:bcc_h3,3)
  real(rt), intent(in   ) :: dx, dy, dz, dt, a_old, a_new

    real(rt)                :: bx, by ,bz, dtdxinv, dtdyinv, dtdzinv, dtinv
  integer                 :: i, j, k

  !------------------------------- Fixing Flux for UEDEN  --------------------------------------------------
  dtdxinv = dt/dx
  dtdyinv = dt/dy
  dtdzinv = dt/dz
  dtinv   = 1.d0/dt
  ! a stuff TODO
  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
!Evolve cell centered non-solenoidal B
        bcc(i,j,k,1:3) = bcc(i,j,k,1:3) + dtdxinv*(flxx(i,j,k,umagx:umagz) - flxx(i+1,j,k,umagx:umagz)) &
                       +  dtdyinv*(flxy(i,j,k,umagx:umagz) - flxy(i,j+1,k,umagx:umagz)) & 
                       +  dtdzinv*(flxz(i,j,k,umagx:umagz) - flxz(i,j,k+1,umagx:umagz))
        bx = 0.5d0*(bxout(i,j,k) + bxout(i+1,j,k)) 
        by = 0.5d0*(byout(i,j,k) + byout(i,j+1,k)) 
        bz = 0.5d0*(bzout(i,j,k) + bzout(i,j,k+1)) 
        hydro_src(i,j,k,ueden) = hydro_src(i,j,k,ueden) + 0.5d0*dtinv*(bx*bx +by*by + bz*bz &
                          - dot_product(bcc(i,j,k,1:3),bcc(i,j,k,1:3)))
       !Corrected Magnetic Energy for solenoidal mag fields
      enddo
    enddo
  enddo
  end subroutine flux_ener_corr

!Combine Fluxes into "hydro_src" made for more general incase we wish to do heating and cooling. 
  subroutine  flux_combo(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                      hydro_src, hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3, &
                      flux1, flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3, &
                      flux2, flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3, &
                      flux3, flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3, &
                      divu_cc, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, lo, hi, &
                      dx, dy, dz, dt, a_old, a_new)                  

       use amrex_fort_module, only : rt => amrex_real 
       use meth_params_module 

  implicit none 
  integer, intent(in    ) :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
  integer, intent(in    ) :: flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3
  integer, intent(in    ) :: flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3
  integer, intent(in    ) :: flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3 
  integer, intent(in    ) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
  integer, intent(in    ) :: hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3
  integer, intent(in    ) :: lo(3), hi(3) 

  real(rt), intent(in   ) :: uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3, NVAR)
  real(rt), intent(in   ) :: flux1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, flux1_l3:flux1_h3, NVAR) 
  real(rt), intent(in   ) :: flux2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, flux2_l3:flux2_h3, NVAR) 
  real(rt), intent(in   ) :: flux3(flux3_l1:flux3_h1, flux3_l2:flux3_h2, flux3_l3:flux3_h3, NVAR) 
  real(rt), intent(inout) :: hydro_src(hsrc_l1:hsrc_h1, hsrc_l2:hsrc_h2, hsrc_l3:hsrc_h3, NVAR) 
  real(rt), intent(in   ) :: divu_cc(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3)
  real(rt), intent(in   ) :: dx, dy, dz, dt, a_old, a_new

  integer :: i, j, k, n 
  real(rt):: dxinv, dyinv, dzinv
  dxinv = 1.d0/dx
  dyinv = 1.d0/dy 
  dzinv = 1.d0/dz

  do n = 1, nvar 
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1) 
            if(n.le.ueden) then 
            hydro_src(i,j,k,n) = (flux1(i,j,k,n) - flux1(i+1,j,k,n))*dxinv &
                               + (flux2(i,j,k,n) - flux2(i,j+1,k,n))*dyinv & 
                               + (flux3(i,j,k,n) - flux3(i,j,k+1,n))*dzinv
            else
            hydro_src(i,j,k,n) = 0.0d0 
!           Species TODO
            endif
!           Source Terms and evolving a TODO 
           enddo
        enddo
     enddo
  enddo         
  end subroutine flux_combo


  subroutine gen_bcc(lo, hi, bcc, bcc_l1, bcc_l2, bcc_l3, bcc_h1, bcc_h2, bcc_h3, &
                   bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
                   byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
                   bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3)
   
       use amrex_fort_module, only : rt => amrex_real 
       use meth_params_module 

  implicit none 
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: bcc_l1 , bcc_l2 , bcc_l3 , bcc_h1 , bcc_h2 , bcc_h3
  integer,  intent(in   ) :: bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
  integer,  intent(in   ) :: byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
  integer,  intent(in   ) :: bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
  
  real(rt), intent(inout) :: bcc (bcc_l1 :bcc_h1 , bcc_l2 :bcc_h2, bcc_l3 :bcc_h3 , 3)
  real(rt), intent(in   ) :: bxin(bxin_l1:bxin_h1,bxin_l2:bxin_h2,bxin_l3:bxin_h3)
  real(rt), intent(in   ) :: byin(byin_l1:byin_h1,byin_l2:byin_h2,byin_l3:byin_h3)
  real(rt), intent(in   ) :: bzin(bzin_l1:bzin_h1,bzin_l2:bzin_h2,bzin_l3:bzin_h3)

  integer                 :: i, j, k 

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          bcc(i,j,k,1) = 0.5d0*(bxin(i,j,k) + bxin(i+1,j,k))
          bcc(i,j,k,2) = 0.5d0*(byin(i,j,k) + byin(i,j+1,k))
          bcc(i,j,k,3) = 0.5d0*(bzin(i,j,k) + bzin(i,j,k+1))
        enddo
     enddo
  enddo

  end subroutine gen_bcc




````

