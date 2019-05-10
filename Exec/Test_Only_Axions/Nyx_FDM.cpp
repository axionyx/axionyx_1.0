#ifdef FDM

#include "Nyx.H"
#include "Nyx_F.H"
#include <AMReX_Particles_F.H>
#include <AMReX_MultiFab.H>
#ifdef GRAVITY
#	include "Gravity.H"
#	include <Gravity_F.H>
#endif

//testing swfft
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <SWFFT_Test_F.H>

// These are for SWFFT
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <string>

#define ALIGN 16

using namespace amrex;


void Nyx::advance_FDM_FD (amrex::Real time,
                      amrex::Real dt,
                      amrex::Real a_old,
                      amrex::Real a_new)
{
    amrex::Real se, ske;
    const amrex::Real prev_time    = state[State_Type].prevTime();
    const amrex::Real cur_time     = state[State_Type].curTime();
    const int  finest_level = parent->finestLevel();

    const amrex::Real a_half = 0.5 * (a_old + a_new);

    amrex::MultiFab&  Ax_old = get_old_data(Axion_Type);
    amrex::MultiFab&  Ax_new = get_new_data(Axion_Type);
    amrex::MultiFab& Phi_old = get_old_data(PhiGrav_Type);
    amrex::MultiFab& Phi_new = get_new_data(PhiGrav_Type);

    if ( amrex::ParallelDescriptor::IOProcessor() ){
	std::cout << "Advancing the axions at level " << level <<  "...\n";
    }
#ifndef NDEBUG
    if (Ax_old.contains_nan(0, Ax_old.nComp(), 0))
    {
        for (int i = 0; i < Ax_old.nComp(); i++)
        {
            if (Ax_old.contains_nan(i,1,0))
            {
                std::cout << "Testing component i for NaNs: " << i << std::endl;
                amrex::Abort("Ax_old has NaNs in this component::advance_FDM_FD()");
            }
        }
    }
#endif

    const amrex::Real* dx      = geom.CellSize();
    amrex::Real        courno  = -1.0e+200;
    //Need the positions of the physical boundaries for the 'sponge' neccessary for apropriate
    //boundary conditions for the Schroedinger-Poisson system (arXiv:gr-qc/0404014v2)
    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* prob_hi = geom.ProbHi();


    // Define the gravity vector so we can pass this to ca_umdrv.
    // MultiFab grav_vector(grids, BL_SPACEDIM, 3, Fab_allocate);
    const auto& dm = get_level(level).get_new_data(State_Type).DistributionMap();
    amrex::MultiFab grav_vector(grids, dm, BL_SPACEDIM, 1);
    grav_vector.setVal(0);


#ifdef GRAVITY
    gravity->get_old_grav_vector(level, grav_vector, time);
    grav_vector.FillBoundary();
    grav_vector.EnforcePeriodicity(geom.periodicity());
    //TODO check if this is ok
    //geom.FillPeriodicBoundary(grav_vector, 0, BL_SPACEDIM);
#endif
    {

      //  // maybe 4 should be NUM_GROW
      // for (FillPatchIterator
      // 	      fpi(*this,  Ax_new, 4, time, Axion_Type,   0, Nyx::NUM_AX),
      // 		       pfpi(*this, Phi_new, 4, time, PhiGrav_Type, 0, 1);
      //       fpi.isValid() && pfpi.isValid();
      //       ++fpi,++pfpi)
      for (amrex::FillPatchIterator
	     fpi(*this,  Ax_new, 2, time, Axion_Type,   0, Nyx::NUM_AX),
	     pfpi(*this, Phi_new, 1, time, PhiGrav_Type, 0, 1);
	     fpi.isValid() && pfpi.isValid();
	     ++fpi,++pfpi)
       {
        // Create FAB for extended grid values (including boundaries) and fill.
        const int  mfi_index = fpi.index();
        const amrex::Box& bx        = fpi.UngrownBox();
        amrex::FArrayBox& axion     = fpi();
        amrex::FArrayBox& axionout  = Ax_new[fpi];
        amrex::FArrayBox& phiold    = pfpi();

        amrex::Real cflLoc = -1.0e+200;

        if (axion.contains_nan())
        {
		amrex::Abort("Nans in state just before fortran call");
	}

	//axionout.copy(axion);
        BL_FORT_PROC_CALL(FORT_ADVANCE_FDM_FD, fort_advance_fdm_fd)
            (&time, bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(axion),
             BL_TO_FORTRAN(axionout),
             BL_TO_FORTRAN(grav_vector[fpi]),
             BL_TO_FORTRAN(phiold),
             dx, prob_lo, prob_hi, &dt,
             &cflLoc, &a_old, &a_new, verbose);
        courno = std::max(courno, cflLoc);
       }
    }

    grav_vector.clear();

    amrex::ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0)
    {
        if (amrex::ParallelDescriptor::IOProcessor())
            std::cout << "OOPS -- EFFECTIVE CFL AT THIS LEVEL " << level
                      << " IS " << courno << '\n';

        amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
    }

#ifndef NDEBUG

    for (int i = 0; i < Ax_new.nComp(); i++)
    {
        if (Ax_new.contains_nan(i, 1, 0))
        {
            std::cout << "Testing component i for NaNs: " << i << std::endl;
            amrex::Abort("Ax_new has NaNs in this component::advance_FDM_FD()");
        }
    }
    int ang = Ax_new.nGrow();

    if ( amrex::ParallelDescriptor::IOProcessor() ){
	std::cout << "Ax_New.nGrow(): " << ang << "\n";
    }

    for (int i = 0; i < Ax_new.nComp(); i++)
    {
        if (Ax_new.contains_nan(i, 1, ang))
        {
            std::cout << "Testing component i for NaNs: " << i << std::endl;
            amrex::Abort("Ax_new has NaNs _in ghostzones_ in this component::advance_FDM_FD()");
        }
    }
#endif


}
//TODO_JENS: The stub I added. This guy should be called in Nyx::advance_FDM when needed (ie. on level 0), prepare the data, and call swfft_solve.
void Nyx::advance_FDM_FFT (amrex::Real time,
			   amrex::Real dt,
			   amrex::Real a_old,
			   amrex::Real a_new)
{


    // *****************************************
    //define constants
    // *****************************************
    Real m_tt = 2.5;
    Real hbaroverm = 0.01917152 / m_tt;
    Real a_half = 0.5 * (a_old + a_new);
    Real pi = 4 * std::atan(1.0);
    Real tpi = 2 * pi;
    Real h = geom.CellSize(0);
    Real hsq = h*h;

    // *****************************************
    // Defining through Axion_Type
    // *****************************************
    MultiFab& Ax_old = get_old_data(Axion_Type);
    MultiFab& Ax_new = get_new_data(Axion_Type);
    const BoxArray& ba = Ax_old.boxArray();
    const DistributionMapping& dm = Ax_old.DistributionMap();
    MultiFab real_old(ba, dm, 1, 0);
    MultiFab imag_old(ba, dm, 1, 0);
    MultiFab::Copy(real_old, Ax_old, Nyx::AxRe, 0, 1, 0);
    MultiFab::Copy(imag_old, Ax_old, Nyx::AxIm, 0, 1, 0);
    MultiFab real_new(ba, dm, 1, 0);
    MultiFab imag_new(ba, dm, 1, 0);

#ifndef NDEBUG
    if (Ax_old.contains_nan(0, Ax_old.nComp(), 0))
    {
        for (int i = 0; i < Ax_old.nComp(); i++)
        {
            if (Ax_old.contains_nan(i,1,0))
            {
                std::cout << "Testing component i for NaNs: " << i << std::endl;
                amrex::Abort("Ax_old has NaNs in this component::advance_FDM_FFT()");
            }
        }
    }
#endif

    // *****************************************
    //Defining through PhiGrav_Type
    // *****************************************
    MultiFab& grav_phi = get_old_data(PhiGrav_Type);
    MultiFab phi(grav_phi.boxArray(), grav_phi.DistributionMap(), 1, 0);
    MultiFab::Copy(phi, grav_phi, 0, 0, 1, 0);

#ifndef NDEBUG
    if (phi.contains_nan(0, phi.nComp(), 0))
    {
        for (int i = 0; i < phi.nComp(); i++)
        {
            if (phi.contains_nan(i,1,0))
            {
                std::cout << "Testing component phi for NaNs: " << i << std::endl;
                amrex::Abort("phi has NaNs in this component::advance_FDM_FFT()");
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
        int local_index = k*nbx*nby + j*nbx + i;

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


////////////////////////////////////////////////////////////////////////////////
    const int *self = dfft.self_kspace();
    const int *local_ng = dfft.local_ng_kspace();
    const int *global_ng = dfft.global_ng();

    for (MFIter mfi(real_old); mfi.isValid(); ++mfi){
         const Box& box = mfi.validbox();

         size_t local_indx = 0;

         double center_phi = 0.512/2;
         // for(size_t i=0; i<(size_t)local_ng[0]; i++) {
         //     double global_x = double(local_ng[0]*self[0] + i)/double(global_ng[0])*0.512;
         //     for(size_t j=0; j<(size_t)local_ng[1]; j++) {
         //         // double global_y = double(local_ng[1]*self[1] + j)/double(global_ng[1]);
         //         for(size_t k=0; k<(size_t)local_ng[2]; k++) {
         //             // double global_z = double(local_ng[2]*self[2] + k)/double(global_ng[1]);
         //
         //             phi[mfi].dataPtr()[local_indx] = 0.5*std::pow(0.1,-4)*std::pow(global_x-center_phi,2);
         //             local_indx++;
         //         }
         //     }
         // }

         for(size_t i=0; i<(size_t)gridsize; i++){
             phi[mfi].dataPtr()[i] = 0.0;
         }
     }
////////////////////////////////////////////////////////////////////////////////



    //// *******************************************
    //// kick - one full time step
    //// *******************************************
    //// FORTRAN KICK
    for (MFIter mfi(real_old); mfi.isValid(); ++mfi){
        const Box& box = mfi.validbox();
        fort_kick(box.loVect(), box.hiVect(), real_old[mfi].dataPtr(), imag_old[mfi].dataPtr(), phi[mfi].dataPtr(), &hbaroverm, &dt);
    }

    // //// C++ KICK
    // for (MFIter mfi(real_old); mfi.isValid(); ++mfi){
    //     const Box& box = mfi.validbox();
    //     const std::complex<double> imagi(0.0,1.0);
    //
    //     for(size_t i=0; i<(size_t)gridsize; i++){
    //         std::complex<double> psi(real_old[mfi].dataPtr()[i],imag_old[mfi].dataPtr()[i]);
    //         psi = std::exp( imagi * phi[mfi].dataPtr()[i] / hbaroverm  * dt ) * psi;
    //         real_old[mfi].dataPtr()[i] = std::real(psi);
    //         imag_old[mfi].dataPtr()[i] = std::imag(psi);
    //     }
    // }


    //Uncomment if you want to plot a MULTIFAB into a file
    // writeMultiFabAsPlotFile("MFAB_plot000", real_old, "psi_re");

    //  *******************************************
    //  drift - one full time step
    //  *******************************************
    for (MFIter mfi(real_old,false); mfi.isValid(); ++mfi){

        std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
        std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;

        a.resize(gridsize);
        b.resize(gridsize);

        dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);

        for(size_t i=0; i<(size_t)gridsize; i++){
            a[i] = complex_t(real_old[mfi].dataPtr()[i],imag_old[mfi].dataPtr()[i]);
        }

        //  *******************************************
        //  Compute the forward transform
        //  *******************************************
        dfft.forward(&a[0]);

// ////////////////////////////////////////////////////////////////////////////////
//       const int *self = dfft.self_kspace();
//       const int *local_ng = dfft.local_ng_kspace();
//       const int *global_ng = dfft.global_ng();
// ////////////////////////////////////////////////////////////////////////////////



        const std::complex<double> imagi(0.0,1.0);
        size_t local_indx = 0;


        // ///// the cos(k)-1 drift loop
        // for(size_t i=0; i<(size_t)local_ng[0]; i++) {
        //     size_t global_i = local_ng[0]*self[0] + i;
        //
        //     for(size_t j=0; j<(size_t)local_ng[1]; j++) {
        //         size_t global_j = local_ng[1]*self[1] + j;
        //
        //         for(size_t k=0; k<(size_t)local_ng[2]; k++) {
        //             size_t global_k = local_ng[2]*self[2] + k;
        //
        //             double laplace_k = 2./hsq * ( (cos(tpi*double(global_i)/double(global_ng[0])) - 1.) +
        //             (cos(tpi*double(global_j)/double(global_ng[1])) - 1.) +
        //             (cos(tpi*double(global_k)/double(global_ng[2])) - 1.) );
        //
        //             a[local_indx] *= std::exp( imagi * hbaroverm * laplace_k / a_half / a_half / 2.0  * dt );
        //             local_indx++;
        //         }
        //     }
        // }



        //////// k^2 drift loop
        for(size_t i=0; i<(size_t)local_ng[0]; i++) {
            int global_i = local_ng[0]*self[0] + i;
            if (global_i >= global_ng[0]/2){
                global_i = global_i - global_ng[0];
            }


            for(size_t j=0; j<(size_t)local_ng[1]; j++) {
                int global_j = local_ng[1]*self[1] + j;
                if (global_j >= global_ng[1]/2){
                    global_j = global_j - global_ng[1];
                }


                for(size_t k=0; k<(size_t)local_ng[2]; k++) {
                    int global_k = local_ng[2]*self[2] + k;
                    if (global_k >= global_ng[2]/2){
                        global_k = global_k - global_ng[2];
                    }

                    double kx = tpi * double(global_i)/double(global_ng[0]);
                    double ky = tpi * double(global_j)/double(global_ng[1]);
                    double kz = tpi * double(global_k)/double(global_ng[2]);
                    double k2 = (kx*kx + ky*ky + kz*kz)/h/h;

                    a[local_indx] *= std::exp(- imagi * hbaroverm * k2 / a_half / a_half / 2.0  * dt );

                    local_indx++;
                }
            }
        }

        // *******************************************
        // Compute the backward transformdfft.global_size
        // *******************************************
        dfft.backward(&a[0]);

        size_t global_size  = dfft.global_size();
        for(size_t i=0; i<(size_t)gridsize; i++) {
            real_new[mfi].dataPtr()[i] = std::real(a[i])/global_size;
            imag_new[mfi].dataPtr()[i] = std::imag(a[i])/global_size;
        }
    }//MFIter loop

    //  *******************************************
    //  Update axion state
    //  *******************************************
    Ax_new.ParallelCopy(real_new, 0, Nyx::AxRe, 1, real_new.nGrow(), Ax_new.nGrow(), geom.periodicity(),FabArrayBase::COPY);
    Ax_new.ParallelCopy(imag_new, 0, Nyx::AxIm, 1, imag_new.nGrow(), Ax_new.nGrow(), geom.periodicity(),FabArrayBase::COPY);

    for (MFIter mfi(Ax_new); mfi.isValid(); ++mfi){
        const Box& box = Ax_new[mfi].box();
        fort_ax_fields(Ax_new[mfi].dataPtr(), box.loVect(), box.hiVect());
    }

#ifndef NDEBUG
    if (Ax_new.contains_nan(0, Ax_new.nComp(), 0))
    {
        for (int i = 0; i < Ax_new.nComp(); i++)
        {
            if (Ax_new.contains_nan(i,1,0))
            {
                std::cout << "Testing component i for NaNs: " << i << std::endl;
                amrex::Abort("Ax_new has NaNs in this component::advance_FDM_FFT()");
            }
        }
    }
#endif

}
#endif //FDM
