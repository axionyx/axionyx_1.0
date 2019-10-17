
# File advance\_FDM\_PS.cpp

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**advance\_FDM\_PS.cpp**](advance__FDM__PS_8cpp.md)

[Go to the documentation of this file.](advance__FDM__PS_8cpp.md) 


````cpp
#ifdef FDM

#include <AMReX_BLProfiler.H>
#include "Nyx.H"
#include "Nyx_F.H"
//#include <AMReX_Particles_F.H>
#include <AMReX_MultiFab.H>
#ifdef GRAVITY
#include "Gravity.H"
#include <Gravity_F.H>
#endif

//testing swfft
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
//#include <SWFFT_Test_F.H>

// These are for SWFFT
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>
#include "Stopwatch.H"

#include <string>

#define ALIGN 16

using namespace amrex;

void copy_c2fortran(std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >& a, 
            MultiFab &refab, MultiFab &imfab, MFIter& mfi);

void copy_fortran2c(std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >& a, 
            MultiFab &refab, MultiFab &imfab, MFIter& mfi);

void drift(hacc::Dfft &dfft, MultiFab &Ax_new, int const gridsize, Real const dt, Real const h, Real const a_half, 
       Real const hbaroverm, std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* a,std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* b);

void fdm_timestep(hacc::Dfft &dfft, MultiFab &Ax_new, MultiFab &phi,  Gravity::Gravity *gravity, Geometry &geom,
          int const level, int const gridsize, Real const h, Real const dt_c,  
          Real const a_c, Real const dt_d, Real const a_d, Real const a_new, 
          Real const hbaroverm,std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* a,std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* b);

// //TODO_JENS: The stub I added. This guy should be called in Nyx::advance_FDM when needed (ie. on level 0), prepare the data, and call swfft_solve.
// void Nyx::advance_FDM_FFT (amrex::Real time,
//             amrex::Real dt,
//             amrex::Real a_old,
//             amrex::Real a_new)
// {


//     // *****************************************
//     //define constants
//     // *****************************************
//     Real m_tt = 2.5;
//     Real hbaroverm = 0.01917152 / m_tt;
//     Real a_half = 0.5 * (a_old + a_new);
//     Real pi = 4 * std::atan(1.0);
//     Real tpi = 2 * pi;
//     Real h = geom.CellSize(0);
//     Real hsq = h*h;

//     // *****************************************
//     // Defining through Axion_Type
//     // *****************************************
//     MultiFab& Ax_old = get_old_data(Axion_Type);
//     MultiFab& Ax_new = get_new_data(Axion_Type);
//     const BoxArray& ba = Ax_old.boxArray();
//     const DistributionMapping& dm = Ax_old.DistributionMap();
//     MultiFab real_old(ba, dm, 1, 0);
//     MultiFab imag_old(ba, dm, 1, 0);
//     MultiFab::Copy(real_old, Ax_old, Nyx::AxRe, 0, 1, 0);
//     MultiFab::Copy(imag_old, Ax_old, Nyx::AxIm, 0, 1, 0);
//     MultiFab real_new(ba, dm, 1, 0);
//     MultiFab imag_new(ba, dm, 1, 0);

// #ifndef NDEBUG
//     if (Ax_old.contains_nan(0, Ax_old.nComp(), 0))
//     {
//         for (int i = 0; i < Ax_old.nComp(); i++)
//         {
//             if (Ax_old.contains_nan(i,1,0))
//             {
//                 std::cout << "Testing component i for NaNs: " << i << std::endl;
//                 amrex::Abort("Ax_old has NaNs in this component::advance_FDM_FFT()");
//             }
//         }
//     }
// #endif

//     // *****************************************
//     //Defining through PhiGrav_Type
//     // *****************************************
//     MultiFab& grav_phi = get_old_data(PhiGrav_Type);
//     MultiFab phi(grav_phi.boxArray(), grav_phi.DistributionMap(), 1, 0);
//     MultiFab::Copy(phi, grav_phi, 0, 0, 1, 0);

// #ifndef NDEBUG
//     if (phi.contains_nan(0, phi.nComp(), 0))
//     {
//         for (int i = 0; i < phi.nComp(); i++)
//         {
//             if (phi.contains_nan(i,1,0))
//             {
//                 std::cout << "Testing component phi for NaNs: " << i << std::endl;
//                 amrex::Abort("phi has NaNs in this component::advance_FDM_FFT()");
//             }
//         }
//     }
// #endif

//     // *****************************************
//     // We assume that all grids have the same size hence
//     // we have the same nx,ny,nz on all ranks
//     // *****************************************
//     int nx = ba[0].size()[0];
//     int ny = ba[0].size()[1];
//     int nz = ba[0].size()[2];
//     int gridsize = nx*ny*nz;

//     Box domain(geom.Domain());

//     int nbx = domain.length(0) / nx;
//     int nby = domain.length(1) / ny;
//     int nbz = domain.length(2) / nz;
//     int nboxes = nbx * nby * nbz;

//     if (nboxes != ba.size())
//         amrex::Error("NBOXES NOT COMPUTED CORRECTLY");

//     // *****************************************
//     // This unfortunately seems neccessary (cf. amrex/Src/Extern/SWFFT/README)
//     // *****************************************
//     if(ParallelDescriptor::NProcs() != nboxes)
//         amrex::Error("Number of MPI ranks has to equal number of root level grids!");

//     Vector<int> rank_mapping;
//     rank_mapping.resize(nboxes);

//     for (int ib = 0; ib < nboxes; ++ib)
//     {
//         int i = ba[ib].smallEnd(0) / nx;
//         int j = ba[ib].smallEnd(1) / ny;
//         int k = ba[ib].smallEnd(2) / nz;
//         int local_index = k*nbx*nby + j*nbx + i;

//         rank_mapping[local_index] = dm[ib];

//         if (verbose)
//             amrex::Print() << "LOADING RANK NUMBER " << dm[ib] << " FOR GRID NUMBER " << ib
//                  << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
//     }

//     // *****************************************
//     // Assume for now that nx = ny = nz
//     // *****************************************
//     int Ndims[3] = { nbz, nby, nbx };
//     int     n[3] = {domain.length(2), domain.length(1), domain.length(0)};
//     hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
//     hacc::Dfft dfft(d);



//     const int *self = dfft.self_kspace();
//     const int *local_ng = dfft.local_ng_kspace();
//     const int *global_ng = dfft.global_ng();

//     // for (MFIter mfi(real_old); mfi.isValid(); ++mfi){
//     //      const Box& box = mfi.validbox();

//     //      size_t local_indx = 0;

//     //      double center_phi = 0.512/2;
//     //      // for(size_t i=0; i<(size_t)local_ng[0]; i++) {
//     //      //     double global_x = double(local_ng[0]*self[0] + i)/double(global_ng[0])*0.512;
//     //      //     for(size_t j=0; j<(size_t)local_ng[1]; j++) {
//     //      //         // double global_y = double(local_ng[1]*self[1] + j)/double(global_ng[1]);
//     //      //         for(size_t k=0; k<(size_t)local_ng[2]; k++) {
//     //      //             // double global_z = double(local_ng[2]*self[2] + k)/double(global_ng[1]);
//     //      //
//     //      //             phi[mfi].dataPtr()[local_indx] = 0.5*std::pow(0.1,-4)*std::pow(global_x-center_phi,2);
//     //      //             local_indx++;
//     //      //         }
//     //      //     }
//     //      // }

//     //      for(size_t i=0; i<(size_t)gridsize; i++){
//     //          phi[mfi].dataPtr()[i] = 0.0;
//     //      }
//     //  }
// ////////////////////////////////////////////////////////////////////////////////


//     //// *******************************************
//     //// kick - one full time step
//     //// *******************************************
//     //// FORTRAN KICK
//     for (MFIter mfi(real_old); mfi.isValid(); ++mfi){
//         const Box& box = mfi.validbox();
//         fort_kick(box.loVect(), box.hiVect(), real_old[mfi].dataPtr(), imag_old[mfi].dataPtr(), phi[mfi].dataPtr(), &hbaroverm, &dt);
//     }

//     // //// C++ KICK
//     // for (MFIter mfi(real_old); mfi.isValid(); ++mfi){
//     //     const Box& box = mfi.validbox();
//     //     const std::complex<double> imagi(0.0,1.0);
//     //
//     //     for(size_t i=0; i<(size_t)gridsize; i++){
//     //         std::complex<double> psi(real_old[mfi].dataPtr()[i],imag_old[mfi].dataPtr()[i]);
//     //         psi = std::exp( imagi * phi[mfi].dataPtr()[i] / hbaroverm  * dt ) * psi;
//     //         real_old[mfi].dataPtr()[i] = std::real(psi);
//     //         imag_old[mfi].dataPtr()[i] = std::imag(psi);
//     //     }
//     // }


//     //Uncomment if you want to plot a MULTIFAB into a file
//     // writeMultiFabAsPlotFile("MFAB_plot000", real_old, "psi_re");

//     //  *******************************************
//     //  drift - one full time step
//     //  *******************************************
//     for (MFIter mfi(real_old,false); mfi.isValid(); ++mfi){

//         std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
//         std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;

//         a.resize(gridsize);
//         b.resize(gridsize);

//         dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);

//         for(size_t i=0; i<(size_t)gridsize; i++){
//             a[i] = complex_t(real_old[mfi].dataPtr()[i],imag_old[mfi].dataPtr()[i]);
//         }

//         //  *******************************************
//         //  Compute the forward transform
//         //  *******************************************
//         dfft.forward(&a[0]);


//         const std::complex<double> imagi(0.0,1.0);
//         size_t local_indx = 0;


//         // ///// the cos(k)-1 drift loop
//         // for(size_t i=0; i<(size_t)local_ng[0]; i++) {
//         //     size_t global_i = local_ng[0]*self[0] + i;
//         //
//         //     for(size_t j=0; j<(size_t)local_ng[1]; j++) {
//         //         size_t global_j = local_ng[1]*self[1] + j;
//         //
//         //         for(size_t k=0; k<(size_t)local_ng[2]; k++) {
//         //             size_t global_k = local_ng[2]*self[2] + k;
//         //
//         //             double laplace_k = 2./hsq * ( (cos(tpi*double(global_i)/double(global_ng[0])) - 1.) +
//         //             (cos(tpi*double(global_j)/double(global_ng[1])) - 1.) +
//         //             (cos(tpi*double(global_k)/double(global_ng[2])) - 1.) );
//         //
//         //             a[local_indx] *= std::exp( imagi * hbaroverm * laplace_k / a_half / a_half / 2.0  * dt );
//         //             local_indx++;
//         //         }
//         //     }
//         // }



//         //////// k^2 drift loop
//         for(size_t i=0; i<(size_t)local_ng[0]; i++) {
//             int global_i = local_ng[0]*self[0] + i;
//             if (global_i >= global_ng[0]/2){
//                 global_i = global_i - global_ng[0];
//             }


//             for(size_t j=0; j<(size_t)local_ng[1]; j++) {
//                 int global_j = local_ng[1]*self[1] + j;
//                 if (global_j >= global_ng[1]/2){
//                     global_j = global_j - global_ng[1];
//                 }


//                 for(size_t k=0; k<(size_t)local_ng[2]; k++) {
//                     int global_k = local_ng[2]*self[2] + k;
//                     if (global_k >= global_ng[2]/2){
//                         global_k = global_k - global_ng[2];
//                     }

//                     double kx = tpi * double(global_i)/double(global_ng[0]);
//                     double ky = tpi * double(global_j)/double(global_ng[1]);
//                     double kz = tpi * double(global_k)/double(global_ng[2]);
//                     double k2 = (kx*kx + ky*ky + kz*kz)/h/h;

//                     a[local_indx] *= std::exp(- imagi * hbaroverm * k2 / a_half / a_half / 2.0  * dt );

//                     local_indx++;
//                 }
//             }
//         }

//         // *******************************************
//         // Compute the backward transformdfft.global_size
//         // *******************************************
//         dfft.backward(&a[0]);

//         size_t global_size  = dfft.global_size();
//         for(size_t i=0; i<(size_t)gridsize; i++) {
//             real_new[mfi].dataPtr()[i] = std::real(a[i])/global_size;
//             imag_new[mfi].dataPtr()[i] = std::imag(a[i])/global_size;
//         }
//     }//MFIter loop

//     //  *******************************************
//     //  Update axion state
//     //  *******************************************
//     Ax_new.ParallelCopy(real_new, 0, Nyx::AxRe, 1, real_new.nGrow(), Ax_new.nGrow(), geom.periodicity(),FabArrayBase::COPY);
//     Ax_new.ParallelCopy(imag_new, 0, Nyx::AxIm, 1, imag_new.nGrow(), Ax_new.nGrow(), geom.periodicity(),FabArrayBase::COPY);

//     for (MFIter mfi(Ax_new); mfi.isValid(); ++mfi){
//         const Box& box = Ax_new[mfi].box();
//         fort_ax_fields(Ax_new[mfi].dataPtr(), box.loVect(), box.hiVect());
//     }

// #ifndef NDEBUG
//     if (Ax_new.contains_nan(0, Ax_new.nComp(), 0))
//     {
//         for (int i = 0; i < Ax_new.nComp(); i++)
//         {
//             if (Ax_new.contains_nan(i,1,0))
//             {
//                 std::cout << "Testing component i for NaNs: " << i << std::endl;
//                 amrex::Abort("Ax_new has NaNs in this component::advance_FDM_FFT()");
//             }
//         }
//     }
// #endif

// }

//TODO_JENS: The stub I added. This guy should be called in Nyx::advance_FDM when needed (ie. on level 0), prepare the data, and call swfft_solve.
void Nyx::advance_FDM_PS(amrex::Real time,
             amrex::Real dt,
             amrex::Real a_old,
             amrex::Real a_new)
{

  BL_PROFILE("Nyx::advance_FDM_PS()");
  Stopwatch::startlap("advance_FDM_PS",1);
  Stopwatch::startlap("init",2);

    // *****************************************
    //define constants
    // *****************************************
    const Real h = geom.CellSize(0);

    Real w0,w1,w2,w3;
    Real c1,c2,c3,c4,c5,c6,c7,c8;
    Real d1,d2,d3,d4,d5,d6,d7,d8;
    Real a_c1,a_c2,a_c3,a_c4,a_c5,a_c6,a_c7,a_c8;
    Real a_d1,a_d2,a_d3,a_d4,a_d5,a_d6,a_d7,a_d8;
    Real a_time;

    //******************************************
    //define weights for time steps (see Levkov et. al. 2018)
    //******************************************

    if(order==6){
      w1 = -1.17767998417887;
      w2 = 0.235573213359359;
      w3 = 0.784513610477560;
      w0 = 1.0-2.0*(w1+w2+w3);

      c1 = w3/2.0*dt;
      c2 = (w2+w3)/2.0*dt;
      c3 = (w1+w2)/2.0*dt;
      c4 = (w0+w1)/2.0*dt;
      c5 = (w0+w1)/2.0*dt;
      c6 = (w1+w2)/2.0*dt;
      c7 = (w2+w3)/2.0*dt;
      c8 = w3/2.0*dt;

      d1 = w3*dt;
      d2 = w2*dt;
      d3 = w1*dt;
      d4 = w0*dt;
      d5 = w1*dt;
      d6 = w2*dt;
      d7 = w3*dt;
      d8 = 0.0;

      a_time = state[Axion_Type].prevTime()+0.5*c1;
      a_c1 = get_comoving_a(a_time);
      a_time += 0.5*(c1+c2);
      a_c2 = get_comoving_a(a_time);
      a_time += 0.5*(c2+c3);
      a_c3 = get_comoving_a(a_time);
      a_time += 0.5*(c3+c4);
      a_c4 = get_comoving_a(a_time);
      a_time += 0.5*(c4+c5);
      a_c5 = get_comoving_a(a_time);
      a_time += 0.5*(c5+c6);
      a_c6 = get_comoving_a(a_time);
      a_time += 0.5*(c6+c7);
      a_c7 = get_comoving_a(a_time);
      a_time += 0.5*(c7+c8);
      a_c8 = get_comoving_a(a_time);

      a_time = state[Axion_Type].prevTime()+0.5*d1;
      a_d1 = get_comoving_a(a_time);
      a_time += 0.5*(d1+d2);
      a_d2 = get_comoving_a(a_time);
      a_time += 0.5*(d2+d3);
      a_d3 = get_comoving_a(a_time);
      a_time += 0.5*(d3+d4);
      a_d4 = get_comoving_a(a_time);
      a_time += 0.5*(d4+d5);
      a_d5 = get_comoving_a(a_time);
      a_time += 0.5*(d5+d6);
      a_d6 = get_comoving_a(a_time);
      a_time += 0.5*(d6+d7);
      a_d7 = get_comoving_a(a_time);
      a_time += 0.5*(d7+d8);
      a_d8 = get_comoving_a(a_time);

    }else if(order==2){

    c1 = 0.5*dt;
    c2 = 0.5*dt;

    d1 = dt;
    d2 = 0.0;

    a_time = state[Axion_Type].prevTime()+0.5*c1;
    a_c1 = get_comoving_a(a_time);
    a_time += 0.5*(c1+c2);
    a_c2 = get_comoving_a(a_time);

    a_time = state[Axion_Type].prevTime()+0.5*d1;
    a_d1 = get_comoving_a(a_time);
    a_time += 0.5*(d1+d2);
    a_d2 = get_comoving_a(a_time);

    }else
      amrex::Error("Order of algorithm not implemented!");

    // *****************************************
    // Defining through Axion_Type
    // *****************************************
    MultiFab& Ax_old = get_old_data(Axion_Type);
    MultiFab& Ax_new = get_new_data(Axion_Type);
    const BoxArray& ba = Ax_old.boxArray();
    const DistributionMapping& dm = Ax_old.DistributionMap();

#ifndef NDEBUG
    if (Ax_old.contains_nan(0, Ax_old.nComp(), Ax_old.nGrow()))
    {
        for (int i = 0; i < Ax_old.nComp(); i++)
        {
      if (Ax_old.contains_nan(i,1,Ax_old.nGrow()))
            {
                std::cout << "Testing component i for NaNs: " << i << std::endl;
                amrex::Abort("Ax_old has NaNs in this component::advance_FDM_FFT_fourth_order()");
            }
        }
    }
#endif

    // *****************************************
    //Defining through PhiGrav_Type
    // *****************************************
    MultiFab& grav_phi = get_old_data(PhiGrav_Type);
    MultiFab phi(grav_phi.boxArray(), grav_phi.DistributionMap(), 1, grav_phi.nGrow());
    MultiFab::Copy(phi, grav_phi, 0, 0, 1, grav_phi.nGrow());

#ifndef NDEBUG
    if (phi.contains_nan(0, phi.nComp(), phi.nGrow()))
    {
        for (int i = 0; i < phi.nComp(); i++)
        {
      if (phi.contains_nan(i,1,phi.nGrow()))
            {
                std::cout << "Testing component phi for NaNs: " << i << std::endl;
                amrex::Abort("phi has NaNs in this component::advance_FDM_FFT_fourth_order()");
            }
        }
    }
#endif

    // *****************************************
    // We assume that all grids have the same size hence
    // we have the same nx,ny,nz on all ranks
    // *****************************************
    int nx = ba[0].size()[0];
    int ny = ba[0].size()[1];
    int nz = ba[0].size()[2];
    int gridsize = nx*ny*nz;

    Box domain(geom.Domain());

    int nbx = domain.length(0) / nx;
    int nby = domain.length(1) / ny;
    int nbz = domain.length(2) / nz;
    int nboxes = nbx * nby * nbz;

    if (nboxes != ba.size())
        amrex::Error("NBOXES NOT COMPUTED CORRECTLY");

    // *****************************************
    // This unfortunately seems neccessary (cf. amrex/Src/Extern/SWFFT/README)
    // *****************************************
    if(ParallelDescriptor::NProcs() != nboxes)
        amrex::Error("Number of MPI ranks has to equal number of root level grids!");

    Vector<int> rank_mapping;
    rank_mapping.resize(nboxes);

    for (int ib = 0; ib < nboxes; ++ib)
    {
        int i = ba[ib].smallEnd(0) / nx;
        int j = ba[ib].smallEnd(1) / ny;
        int k = ba[ib].smallEnd(2) / nz;
        int local_index = i*nbx*nby + j*nbx + k;

        rank_mapping[local_index] = dm[ib];

        if (verbose)
            amrex::Print() << "LOADING RANK NUMBER " << dm[ib] << " FOR GRID NUMBER " << ib
                 << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
    }
    Stopwatch::stoplap();
    Stopwatch::startlap("dfft prep",2);
    // *****************************************
    // Assume for now that nx = ny = nz
    // *****************************************
    int Ndims[3] = { nbz, nby, nbx };
    int     n[3] = {domain.length(2), domain.length(1), domain.length(0)};
    hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
    hacc::Dfft dfft(d);

    //  *******************************************
    //  prepare fft
    //  *******************************************

    std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
    std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;
    a.resize(gridsize);
    b.resize(gridsize);
    dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);

    for (MFIter mfi(Ax_old,false); mfi.isValid(); ++mfi){
      const Array4<Real> const& arr = Ax_old.array(mfi);
      const Box& bx = mfi.validbox();
      const Dim3 lo = amrex::lbound(bx);
      const Dim3 hi = amrex::ubound(bx);
      const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};
#ifdef _OPENMP
#pragma omp parallel for       
#endif
      for(size_t i=0; i<(size_t)w.x; i++) {
    for(size_t j=0; j<(size_t)w.y; j++) {
      for(size_t k=0; k<(size_t)w.z; k++) {
        size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*i+(size_t)w.z*j+k;
        complex_t temp(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxRe),arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxIm));
        a[local_indx_threaded] = temp;
      }}}
    }

    Stopwatch::stoplap();
    //  *******************************************
    //  higher order time steps
    //  *******************************************
  Stopwatch::startlap("steps",2);
  if(order==6){
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c1, a_c1, d1, a_d1, a_new, hbaroverm, &a, &b);
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c2, a_c2, d2, a_d2, a_new, hbaroverm, &a, &b);
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c3, a_c3, d3, a_d3, a_new, hbaroverm, &a, &b);
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c4, a_c4, d4, a_d4, a_new, hbaroverm, &a, &b);
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c5, a_c5, d5, a_d5, a_new, hbaroverm, &a, &b);
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c6, a_c6, d6, a_d6, a_new, hbaroverm, &a, &b);
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c7, a_c7, d7, a_d7, a_new, hbaroverm, &a, &b);
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c8, a_c8, d8, a_d8, a_new, hbaroverm, &a, &b);
  }else if(order==2){
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c1, a_c1, d1, a_d1, a_new, hbaroverm, &a, &b);
    fdm_timestep(dfft, Ax_new, phi, gravity, geom, level, gridsize, h, c2, a_c2, d2, a_d2, a_new, hbaroverm, &a, &b);
  }

  Stopwatch::stoplap("steps");

    //  *******************************************
    //  copy everything back to Axion_State
    //  *******************************************
    Stopwatch::startlap("finalize",2);

  for (MFIter mfi(Ax_new,false); mfi.isValid(); ++mfi){
    Array4<Real> const& arr = Ax_new[mfi].array();
    const Box& bx = mfi.validbox();
    const Dim3 lo = amrex::lbound(bx);
    const Dim3 hi = amrex::ubound(bx);
    const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};
#ifdef _OPENMP
#pragma omp parallel for       
#endif
    for(size_t i=0; i<(size_t)w.x; i++) {
      for(size_t j=0; j<(size_t)w.y; j++) {
    for(size_t k=0; k<(size_t)w.z; k++) {
      size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*i+(size_t)w.z*j+k;
      complex_t temp = a[local_indx_threaded];
      arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxRe)=std::real(temp);
      arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxIm)=std::imag(temp);
      arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens)=std::real(temp)*std::real(temp)+std::imag(temp)*std::imag(temp);
        }}}
  }

  Ax_new.FillBoundary(geom.periodicity());

    
#ifndef NDEBUG
  if (Ax_new.contains_nan(0, Ax_new.nComp(), Ax_new.nGrow()))
      {
        for (int i = 0; i < Ax_new.nComp(); i++)
      {
            if (Ax_new.contains_nan(i,1,Ax_new.nGrow()))
          {
                std::cout << "Testing component i for NaNs: " << i << std::endl;
                amrex::Abort("Ax_new has NaNs in this component::advance_FDM_FFT_fourth_order()");
          }
      }
      }
#endif
Stopwatch::stoplap();
Stopwatch::stoplap("advance_FDM_PS");
}

inline void fdm_timestep(hacc::Dfft &dfft, MultiFab &Ax_new, MultiFab &phi,  Gravity::Gravity *gravity, Geometry &geom,
             int const level, int const gridsize, Real const h, Real const dt_c,  
             Real const a_c, Real const dt_d, Real const a_d, Real const a_new, 
             Real const hbaroverm,std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* a,std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* b)
{
  //  *******************************************
  //  drift by dt_c
  //  *******************************************
  Stopwatch::startlap("fdm_timestep",2);
  Stopwatch::startlap("drift",3);
  drift(dfft, Ax_new, gridsize, dt_c, h, a_c, hbaroverm, a, b);
  Ax_new.FillBoundary(geom.periodicity());

#ifndef NDEBUG
  if (Ax_new.contains_nan(Nyx::AxDens,1,Ax_new.nGrow()))
    {
      std::cout << "Testing component i for NaNs: " << Nyx::AxDens << std::endl;
      amrex::Abort("Ax_new has NaNs in this component::fdm_timestep()");
    }
#endif

  Stopwatch::stoplap();
  if(!dt_d) 
  {
      Stopwatch::stoplap("fdm_timestep");
      return;
  }

  //  *******************************************
  //  re-calculate potential
  //  *******************************************

  Stopwatch::startlap("potential",3);
  int fill_interior = 0;
  int grav_n_grow = 0;
  gravity->solve_for_new_phi(level,phi,
                 gravity->get_grad_phi_curr(level),
                 fill_interior, grav_n_grow);
  Stopwatch::stoplap();
  //  *******************************************
  //  kick by dt_d 
  //  *******************************************
  Stopwatch::startlap("kick",3);
  const std::complex<double> imagi(0.0,1.0);
  
  for (MFIter mfi(phi,false); mfi.isValid(); ++mfi){
    Array4<Real> const& arr = phi[mfi].array();
    const Box& bx = mfi.validbox();
    const Dim3 lo = amrex::lbound(bx);
    const Dim3 hi = amrex::ubound(bx);
    const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};
#ifdef _OPENMP
#pragma omp parallel for       
#endif
    for(size_t i=0; i<(size_t)w.x; i++) {
      for(size_t j=0; j<(size_t)w.y; j++) {
        AMREX_PRAGMA_SIMD     
    for(size_t k=0; k<(size_t)w.z; k++) {
      size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*i+(size_t)w.z*j+k;
      (*a)[local_indx_threaded] = std::exp( imagi * arr(i+lo.x,j+lo.y,k+lo.z) * a_new / a_d / hbaroverm  * dt_d ) * (*a)[local_indx_threaded];
        }}}
  }    

  Stopwatch::stoplap();
  Stopwatch::stoplap("fdm_timestep");
}

inline void drift(hacc::Dfft &dfft, MultiFab &Ax_new, 
          int const gridsize, Real const dt, Real const h, Real const a_half, 
          Real const hbaroverm,std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* a,std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* b)
{
  if(!dt) return;

  const Real pi = 4 * std::atan(1.0);
  const Real tpi = 2 * pi;
  const Real hsq = h*h;
  const std::complex<double> imagi(0.0,1.0);
  const int *self = dfft.self_kspace();
  const int *local_ng = dfft.local_ng_kspace();
  const int *global_ng = dfft.global_ng();
  size_t local_indx = 0;

    Stopwatch::startlap("fft",4);
    dfft.forward(&((*a)[0]));    
    Stopwatch::stoplap();
    Stopwatch::startlap("k space operation",4);
#ifdef _OPENMP
#pragma omp parallel for       
#endif
    for(size_t i=0; i<(size_t)local_ng[0]; i++) {
      int global_i = local_ng[0]*self[0] + i;
      if (global_i > global_ng[0]/2.){
    global_i = global_i - global_ng[0];
      }
    
    
      for(size_t j=0; j<(size_t)local_ng[1]; j++) {
    int global_j = local_ng[1]*self[1] + j;
    if (global_j > global_ng[1]/2.){
      global_j = global_j - global_ng[1];
    }
      
        AMREX_PRAGMA_SIMD     
    for(size_t k=0; k<(size_t)local_ng[2]; k++) {
      int global_k = local_ng[2]*self[2] + k;
          size_t local_indx_threaded = (size_t)local_ng[1]*(size_t)local_ng[2]*i+(size_t)local_ng[2]*j+k;
      if (global_k > global_ng[2]/2.){
        global_k = global_k - global_ng[2];
      }
        
      double kx = tpi * double(global_i)/double(global_ng[0]);
      double ky = tpi * double(global_j)/double(global_ng[1]);
      double kz = tpi * double(global_k)/double(global_ng[2]);
      double k2 = (kx*kx + ky*ky + kz*kz)/h/h;
    
      (*a)[local_indx_threaded] *= std::exp(- imagi * hbaroverm * k2 / a_half / a_half / 2.0  * dt );
          (*a)[local_indx_threaded] /= dfft.global_size();  
    }
      }
    }
  Stopwatch::stoplap();
    
  Stopwatch::startlap("fft backwards",4);
    dfft.backward(&((*a)[0]));
  
  Stopwatch::stoplap();
  Stopwatch::startlap("finalize",4);

  for (MFIter mfi(Ax_new,false); mfi.isValid(); ++mfi){
    Array4<Real> const& arr = Ax_new[mfi].array();
    const Box& bx = mfi.validbox();
    const Dim3 lo = amrex::lbound(bx);
    const Dim3 hi = amrex::ubound(bx);
    const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};
#ifdef _OPENMP
#pragma omp parallel for       
#endif
    for(size_t i=0; i<(size_t)w.x; i++) {
      for(size_t j=0; j<(size_t)w.y; j++) {
        AMREX_PRAGMA_SIMD     
    for(size_t k=0; k<(size_t)w.z; k++) {
      size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*i+(size_t)w.z*j+k;
      complex_t temp = (*a)[local_indx_threaded];
      arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens)=std::real(temp)*real(temp)+std::imag(temp)*std::imag(temp);
        }}}

  }
  Stopwatch::stoplap("finalize");
}

inline void copy_fortran2c(std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >& a, 
               MultiFab &refab, MultiFab &imfab, MFIter& mfi)
{
 const Array4<Real> const& imarr = imfab[mfi].array();
 const Array4<Real> const& rearr = refab[mfi].array();
 const Box& bx = mfi.validbox();
 const Dim3 lo = amrex::lbound(bx);
 const Dim3 hi = amrex::ubound(bx);
 const Dim3 w ={hi.x-lo.x,hi.y-lo.y,hi.z-lo.z};

 size_t local_indx = 0;
       for(size_t i=0; i<=(size_t)w.x; i++) {
        for(size_t j=0; j<=(size_t)w.y; j++) {
         for(size_t k=0; k<=(size_t)w.z; k++) {
        complex_t temp(rearr(i+lo.x,j+lo.y,k+lo.z),imarr(i+lo.x,j+lo.y,k+lo.z));
        a[local_indx] = temp;
        local_indx++;
        }}}
}

inline void copy_c2fortran(std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >& a, 
               MultiFab &refab, MultiFab &imfab, MFIter& mfi)
{
 Array4<Real> const& imarr = imfab[mfi].array();
 Array4<Real> const& rearr = refab[mfi].array();
 const Box& bx = mfi.validbox();
 const Dim3 lo = amrex::lbound(bx);
 const Dim3 hi = amrex::ubound(bx);
 const Dim3 w ={hi.x-lo.x,hi.y-lo.y,hi.z-lo.z};
 size_t local_indx = 0;
       for(size_t i=0; i<=(size_t)w.x; i++) {
        for(size_t j=0; j<=(size_t)w.y; j++) {
         for(size_t k=0; k<=(size_t)w.z; k++) {
        complex_t temp = a[local_indx];
        rearr(i+lo.x,j+lo.y,k+lo.z)=std::real(temp);
        imarr(i+lo.x,j+lo.y,k+lo.z)=std::imag(temp);
        local_indx++;
        }}}
}

#endif //FDM
````

