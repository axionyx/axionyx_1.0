#ifdef FDM

#include "Nyx.H"
#include "Nyx_F.H"
//#include <AMReX_Particles_F.H>
#include <AMReX_MultiFab.H>
#ifdef GRAVITY
#	include "Gravity.H"
#	include <Gravity_F.H>
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

#include <string>

#define ALIGN 16

using namespace amrex;

void copy_c2fortran(std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >& a, 
		    MultiFab &refab, MultiFab &imfab, MFIter& mfi);

void copy_fortran2c(std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >& a, 
		    MultiFab &refab, MultiFab &imfab, MFIter& mfi);

void drift(hacc::Dfft &dfft, MultiFab &real, MultiFab &imag, MultiFab &dens, 
	   int const gridsize, Real const dt, Real const h, Real const a_half, 
	   Real const hbaroverm);

void fdm_timestep(hacc::Dfft &dfft, MultiFab &real, MultiFab &imag, MultiFab &dens, 
		  MultiFab &Ax_new, MultiFab &phi,  Gravity::Gravity *gravity, Geometry &geom,
		  int const level, int const gridsize, Real const h, Real const dt_c,  
		  Real const a_c, Real const dt_d, Real const a_d, Real const a_new, 
		  Real const hbaroverm);

// //TODO_JENS: The stub I added. This guy should be called in Nyx::advance_FDM when needed (ie. on level 0), prepare the data, and call swfft_solve.
// void Nyx::advance_FDM_FFT (amrex::Real time,
// 			   amrex::Real dt,
// 			   amrex::Real a_old,
// 			   amrex::Real a_new)
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
    // *****************************************
    //define constants
    // *****************************************
    const Real h = geom.CellSize(0);

    //******************************************
    //define weights for time steps (see Levkov et. al. 2018)
    //******************************************
    const Real w1 = -1.17767998417887;
    const Real w2 = 0.235573213359359;
    const Real w3 = 0.784513610477560;
    const Real w0 = 1.0-2.0*(w1+w2+w3);

    const Real c1 = w3/2.0*dt;
    const Real c2 = (w2+w3)/2.0*dt;
    const Real c3 = (w1+w2)/2.0*dt;
    const Real c4 = (w0+w1)/2.0*dt;
    const Real c5 = (w0+w1)/2.0*dt;
    const Real c6 = (w1+w2)/2.0*dt;
    const Real c7 = (w2+w3)/2.0*dt;
    const Real c8 = w3/2.0*dt;

    const Real d1 = w3*dt;
    const Real d2 = w2*dt;
    const Real d3 = w1*dt;
    const Real d4 = w0*dt;
    const Real d5 = w1*dt;
    const Real d6 = w2*dt;
    const Real d7 = w3*dt;
    const Real d8 = 0.0;

    Real a_time = state[Axion_Type].prevTime()+0.5*c1;
    const Real a_c1 = get_comoving_a(a_time);
    a_time += 0.5*(c1+c2);
    const Real a_c2 = get_comoving_a(a_time);
    a_time += 0.5*(c2+c3);
    const Real a_c3 = get_comoving_a(a_time);
    a_time += 0.5*(c3+c4);
    const Real a_c4 = get_comoving_a(a_time);
    a_time += 0.5*(c4+c5);
    const Real a_c5 = get_comoving_a(a_time);
    a_time += 0.5*(c5+c6);
    const Real a_c6 = get_comoving_a(a_time);
    a_time += 0.5*(c6+c7);
    const Real a_c7 = get_comoving_a(a_time);
    a_time += 0.5*(c7+c8);
    const Real a_c8 = get_comoving_a(a_time);

    a_time = state[Axion_Type].prevTime()+0.5*d1;
    const Real a_d1 = get_comoving_a(a_time);
    a_time += 0.5*(d1+d2);
    const Real a_d2 = get_comoving_a(a_time);
    a_time += 0.5*(d2+d3);
    const Real a_d3 = get_comoving_a(a_time);
    a_time += 0.5*(d3+d4);
    const Real a_d4 = get_comoving_a(a_time);
    a_time += 0.5*(d4+d5);
    const Real a_d5 = get_comoving_a(a_time);
    a_time += 0.5*(d5+d6);
    const Real a_d6 = get_comoving_a(a_time);
    a_time += 0.5*(d6+d7);
    const Real a_d7 = get_comoving_a(a_time);
    a_time += 0.5*(d7+d8);
    const Real a_d8 = get_comoving_a(a_time);

    // *****************************************
    // Defining through Axion_Type
    // *****************************************
    MultiFab& Ax_old = get_old_data(Axion_Type);
    MultiFab& Ax_new = get_new_data(Axion_Type);
    const BoxArray& ba = Ax_old.boxArray();
    const DistributionMapping& dm = Ax_old.DistributionMap();
    MultiFab real(ba, dm, 1, 0);
    MultiFab imag(ba, dm, 1, 0);
    MultiFab dens(ba, dm, 1, 0);
    MultiFab::Copy(real, Ax_old, Nyx::AxRe, 0, 1, 0);
    MultiFab::Copy(imag, Ax_old, Nyx::AxIm, 0, 1, 0);
    dens.setVal(0.0);

#ifndef NDEBUG
    if (Ax_old.contains_nan(0, Ax_old.nComp(), 0))
    {
        for (int i = 0; i < Ax_old.nComp(); i++)
        {
            if (Ax_old.contains_nan(i,1,0))
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
    if (phi.contains_nan(0, phi.nComp(), 0))
    {
        for (int i = 0; i < phi.nComp(); i++)
        {
            if (phi.contains_nan(i,1,0))
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

    // *****************************************
    // Assume for now that nx = ny = nz
    // *****************************************
    int Ndims[3] = { nbz, nby, nbx };
    int     n[3] = {domain.length(2), domain.length(1), domain.length(0)};
    hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
    hacc::Dfft dfft(d);

    //  *******************************************
    //  higher order time steps
    //  *******************************************
    fdm_timestep(dfft, real, imag, dens, Ax_new, phi, gravity, geom, level, gridsize, h, c1, a_c1, d1, a_d1, a_new, hbaroverm);
    fdm_timestep(dfft, real, imag, dens, Ax_new, phi, gravity, geom, level, gridsize, h, c2, a_c2, d2, a_d2, a_new, hbaroverm);
    fdm_timestep(dfft, real, imag, dens, Ax_new, phi, gravity, geom, level, gridsize, h, c3, a_c3, d3, a_d3, a_new, hbaroverm);
    fdm_timestep(dfft, real, imag, dens, Ax_new, phi, gravity, geom, level, gridsize, h, c4, a_c4, d4, a_d4, a_new, hbaroverm);
    fdm_timestep(dfft, real, imag, dens, Ax_new, phi, gravity, geom, level, gridsize, h, c5, a_c5, d5, a_d5, a_new, hbaroverm);
    fdm_timestep(dfft, real, imag, dens, Ax_new, phi, gravity, geom, level, gridsize, h, c6, a_c6, d6, a_d6, a_new, hbaroverm);
    fdm_timestep(dfft, real, imag, dens, Ax_new, phi, gravity, geom, level, gridsize, h, c7, a_c7, d7, a_d7, a_new, hbaroverm);
    fdm_timestep(dfft, real, imag, dens, Ax_new, phi, gravity, geom, level, gridsize, h, c8, a_c8, d8, a_d8, a_new, hbaroverm);

    //  *******************************************
    //  copy everything back to Axion_State
    //  *******************************************
    MultiFab::Copy(Ax_new, real, 0, Nyx::AxRe, 1, 0);
    MultiFab::Copy(Ax_new, imag, 0, Nyx::AxIm, 1, 0);
    MultiFab::Copy(Ax_new, dens, 0, Nyx::AxDens, 1, 0);
    Ax_new.FillBoundary(geom.periodicity());
    
#ifndef NDEBUG
    if (Ax_new.contains_nan(0, Ax_new.nComp(), 0))
      {
        for (int i = 0; i < Ax_new.nComp(); i++)
	  {
            if (Ax_new.contains_nan(i,1,0))
	      {
                std::cout << "Testing component i for NaNs: " << i << std::endl;
                amrex::Abort("Ax_new has NaNs in this component::advance_FDM_FFT_fourth_order()");
	      }
	  }
      }
#endif
    
}

inline void fdm_timestep(hacc::Dfft &dfft, MultiFab &real, MultiFab &imag, MultiFab &dens, 
			 MultiFab &Ax_new, MultiFab &phi,  Gravity::Gravity *gravity, Geometry &geom,
			 int const level, int const gridsize, Real const h, Real const dt_c,  
			 Real const a_c, Real const dt_d, Real const a_d, Real const a_new, 
			 Real const hbaroverm)
{
  //  *******************************************
  //  drift by dt_c
  //  *******************************************
  drift(dfft, real, imag, dens, gridsize, dt_c, h, a_c, hbaroverm);
  if(!dt_d) return;

  //  *******************************************
  //  re-calculate potential
  //  *******************************************
  Ax_new.ParallelCopy(dens, 0, Nyx::AxDens, 1, dens.nGrow(), 
		      Ax_new.nGrow(), geom.periodicity(),FabArrayBase::COPY);
  int fill_interior = 0;
  int grav_n_grow = 0;
  gravity->solve_for_new_phi(level,phi,
			     gravity->get_grad_phi_curr(level),
			     fill_interior, grav_n_grow);
  
  //  *******************************************
  //  kick by dt_d 
  //  *******************************************
  MultiFab phi_temp(phi.boxArray(), phi.DistributionMap(), 1, 0);
  MultiFab::Copy(phi_temp, phi, 0, 0, 1, 0);
  const std::complex<double> imagi(0.0,1.0);
  
  for (MFIter mfi(dens,false); mfi.isValid(); ++mfi){
    std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
    a.resize(gridsize);
    for(size_t i=0; i<(size_t)gridsize; i++){
      a[i] = complex_t(real[mfi].dataPtr()[i],imag[mfi].dataPtr()[i]);
      a[i] = std::exp( imagi * phi_temp[mfi].dataPtr()[i] * a_new / a_d / hbaroverm  * dt_d ) * a[i];
      real[mfi].dataPtr()[i] = std::real(a[i]);
      imag[mfi].dataPtr()[i] = std::imag(a[i]);
      dens[mfi].dataPtr()[i] = std::real(a[i])*std::real(a[i])+std::imag(a[i])*std::imag(a[i]);
    }
  }
}

inline void drift(hacc::Dfft &dfft, MultiFab &real, MultiFab &imag, MultiFab &dens, 
		  int const gridsize, Real const dt, Real const h, Real const a_half, 
		  Real const hbaroverm)
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

  for (MFIter mfi(dens,false); mfi.isValid(); ++mfi){
    
    std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
    std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;
    a.resize(gridsize);
    b.resize(gridsize);
    
    dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);
    
   copy_fortran2c(a,real,imag,mfi);
   // for(size_t i=0; i<(size_t)gridsize; i++){
   //   a[i] = complex_t(real[mfi].dataPtr()[i],imag[mfi].dataPtr()[i]);
   // }

    dfft.forward(&a[0]);
      
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
	  
	  
	for(size_t k=0; k<(size_t)local_ng[2]; k++) {
	  int global_k = local_ng[2]*self[2] + k;
	  if (global_k > global_ng[2]/2.){
	    global_k = global_k - global_ng[2];
	  }
	    
	  double kx = tpi * double(global_i)/double(global_ng[0]);
	  double ky = tpi * double(global_j)/double(global_ng[1]);
	  double kz = tpi * double(global_k)/double(global_ng[2]);
	  double k2 = (kx*kx + ky*ky + kz*kz)/h/h;
	
	  a[local_indx] *= std::exp(- imagi * hbaroverm * k2 / a_half / a_half / 2.0  * dt );
          a[local_indx] /= dfft.global_size();	
	  local_indx++;
	}
      }
    }
    
    dfft.backward(&a[0]);
  
    copy_c2fortran(a,real,imag,mfi);
    //for(size_t i=0; i<(size_t)gridsize; i++) {
    //          real[mfi].dataPtr()[i] = std::real(a[i]);
    //          imag[mfi].dataPtr()[i] = std::imag(a[i]);
    //}
    for(size_t i=0; i<(size_t)gridsize; i++) {
      dens[mfi].dataPtr()[i] = real[mfi].dataPtr()[i]*real[mfi].dataPtr()[i]+imag[mfi].dataPtr()[i]*imag[mfi].dataPtr()[i];
    }
  }  
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
