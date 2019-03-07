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
    // Defining through Axion_Type
    // *****************************************
    MultiFab& Ax_old = get_old_data(Axion_Type);
    MultiFab& Ax_new = get_new_data(Axion_Type);

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
    
    const BoxArray& ba = Ax_old.boxArray();
    const DistributionMapping& dm = Ax_old.DistributionMap();
    MultiFab real(ba, dm, 1, 0);
    MultiFab imag(ba, dm, 1, 0);
    MultiFab::Copy(real, Ax_old, Nyx::AxRe, 0, 1, 0);
    MultiFab::Copy(imag, Ax_old, Nyx::AxIm, 0, 1, 0);

    //Uncomment if you want to plot rhs before swfft
    // writeMultiFabAsPlotFile("MFAB_plot000", real, "rhs");

    // Define pi and (two pi) here                                                                                                                                                                             
    Real  pi = 4 * std::atan(1.0);
    Real tpi = 2 * pi;

    // We assume that all grids have the same size hence                                                                                                                                                        
    // we have the same nx,ny,nz on all ranks                                                                                                                                                                  
    int nx = ba[0].size()[0];
    int ny = ba[0].size()[1];
    int nz = ba[0].size()[2];

    Box domain(geom.Domain());

    int nbx = domain.length(0) / nx;
    int nby = domain.length(1) / ny;
    int nbz = domain.length(2) / nz;
    int nboxes = nbx * nby * nbz;
    if (nboxes != ba.size())
      amrex::Error("NBOXES NOT COMPUTED CORRECTLY");
    amrex::Print() << "Number of boxes:\t" << nboxes << std::endl;

    Vector<int> rank_mapping;
    rank_mapping.resize(nboxes);

    for (int ib = 0; ib < nboxes; ++ib)
      {
        int i = ba[ib].smallEnd(0) / nx;
        int j = ba[ib].smallEnd(1) / ny;
        int k = ba[ib].smallEnd(2) / nz;

        // This would be the "correct" local index if the data wasn't being transformed                                                                                                                         
        // int local_index = k*nbx*nby + j*nbx + i;                                                                                                                                                             

        // This is what we pass to dfft to compensate for the Fortran ordering                                                                                                                                
        //      of amrex data in MultiFabs.                                                                                                                                                                    
        int local_index = i*nby*nbz + j*nbz + k;

        rank_mapping[local_index] = dm[ib];
        if (verbose)
	  amrex::Print() << "LOADING RANK NUMBER " << dm[ib] << " FOR GRID NUMBER " << ib
                         << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
      }

    Real h = geom.CellSize(0);
    Real hsq = h*h;

    Real start_time = amrex::second();

    // Assume for now that nx = ny = nz                                                                                                                                                                       
    int Ndims[3] = { nbz, nby, nbx };
    int     n[3] = {domain.length(2), domain.length(1), domain.length(0)};
    hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
    hacc::Dfft dfft(d);

    //some constants which need to be passed from the init file                                                                                                                                              
    double hbar = 1.;
    double k0 = 1.;
    double m = 1.;
    double hbaroverm = hbar/m;

    //get the potential                                                                                                                                                                                    
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

    //loop for the drift                                                                                                                                                                                           
    for (MFIter mfi(real,false); mfi.isValid(); ++mfi)
      {
	int gid = mfi.index();
        size_t local_size  = dfft.local_size();

	std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
	std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;

        a.resize(nx*ny*nz);
        b.resize(nx*ny*nz);

        dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);

        // *******************************************                                                                                                                                                          
        // Copy real data from Rhs into real part of a -- no ghost cells and                                                                                                                                    
        // put into C++ ordering (not Fortran)                                                                                                                                                                  
        // *******************************************
        complex_t zero(0.0, 0.0);
        size_t local_indx = 0;
        for(size_t k=0; k<(size_t)nz; k++) {
	  for(size_t j=0; j<(size_t)ny; j++) {
	    for(size_t i=0; i<(size_t)nx; i++) {

	      complex_t temp(real[mfi].dataPtr()[local_indx],imag[mfi].dataPtr()[local_indx]);
	      a[local_indx] = temp;
	      local_indx++;

	    }
	  }
        }
        // *******************************************                                                                                                                                                        
        // Same thing as above, for the potentail -only need the real part                                                                                                                                     
        // put into C++ ordering (not Fortran)                                                                                                                                                                 
        // *******************************************                                                                                                                                                        
        local_indx = 0;
        for(size_t k=0; k<(size_t)nz; k++) {
	  for(size_t j=0; j<(size_t)ny; j++) {
	    for(size_t i=0; i<(size_t)nx; i++) {

	      complex_t temp(phi[mfi].dataPtr()[local_indx],0.);
	      b[local_indx] = temp;
	      local_indx++;

	    }
	  }
        }



        //  *******************************************                                                                                                                                                       
        //  Compute the forward transform                                                                                                                                                                     
        //  *******************************************
        dfft.forward(&a[0]);

        //  *******************************************                                                                                                                                                       
        //  drift - one full time step                                                                                                                                                                         
        //  *******************************************                                                                                                                                                          
        const int *self = dfft.self_kspace();
        const int *local_ng = dfft.local_ng_kspace();
        const int *global_ng = dfft.global_ng();
        const std::complex<double> imagi(0.0,1.0);
        local_indx = 0;
        for(size_t k=0; k<(size_t)local_ng[0]; k++) {
	  size_t global_k = local_ng[2]*self[0] + k; //maybe need another condition for if (k < local_ng[2]/2)?                                                                                                  

	  for(size_t j=0; j<(size_t)local_ng[1]; j++) {
	    size_t global_j = local_ng[1]*self[1] + j;

	    for(size_t i=0; i<(size_t)local_ng[2]; i++) {
	      size_t global_i = local_ng[0]*self[2] + i;

	      if (global_i == 0 && global_j == 0 && global_k == 0) {
		a[local_indx] = 0;
	      }
	      else {
		double k2 = k0 * k0 * ( double(global_i) * double(global_i)
					+ double(global_j) * double(global_j)
					+ double(global_k) * double(global_k) );

		a[local_indx] *= std::exp( -imagi * hbar * k2/2.0/m  * dt );

	      }
	      local_indx++;
	    }
	  }
        }

        // *******************************************                                                                                                                                                           
        // Compute the backward transformdfft.global_size                                                                                                                                                        
        // *******************************************                                                                                                                                                           
        dfft.backward(&a[0]);
 
        size_t global_size  = dfft.global_size();
	std::cout << "GOBAL SIZE " << global_size << std::endl;
        double fac = hsq / global_size;
        local_indx = 0;
        for(size_t k=0; k<(size_t)nz; k++) {
	  for(size_t j=0; j<(size_t)ny; j++) {
	    for(size_t i=0; i<(size_t)nx; i++) {

	      real[mfi].dataPtr()[local_indx] = std::real(a[local_indx])/global_size;
	      imag[mfi].dataPtr()[local_indx] = std::imag(a[local_indx])/global_size;
	      local_indx++;

	    }
	  }
        }

      }//MFIter loop


    //  *******************************************                                                                                                                                            
    //  kick - one full time step                                                                                                                                                                 
    //  *******************************************                                                                                            
     for (MFIter mfi(real); mfi.isValid(); ++mfi)
      {
	const Box& box = mfi.validbox();

	Real* real_p = real[mfi].dataPtr();
	Real* imag_p = imag[mfi].dataPtr();
	Real* phi_p  = phi[mfi].dataPtr();

	fort_kick(box.loVect(), box.hiVect(), real_p, imag_p, phi_p, &hbaroverm, &dt);
      }

#ifndef NDEBUG
    if (real.contains_nan(0, real.nComp(), 0))
    {
        for (int i = 0; i < real.nComp(); i++)
        {
            if (real.contains_nan(i,1,0))
            {
                std::cout << "Testing component real for NaNs: " << i << std::endl;
                amrex::Abort("real has NaNs in this component::advance_FDM_FFT()");
            }
        }
    }
#endif
#ifndef NDEBUG
    if (imag.contains_nan(0, imag.nComp(), 0))
    {
        for (int i = 0; i < imag.nComp(); i++)
        {
            if (imag.contains_nan(i,1,0))
            {
                std::cout << "Testing component imag for NaNs: " << i << std::endl;
                amrex::Abort("imag has NaNs in this component::advance_FDM_FFT()");
            }
        }
    }
#endif

    //Uncomment if you want to plot lhs after swfft
    // writeMultiFabAsPlotFile("MFAB_plot111", real, "rhs");

    Ax_new.ParallelCopy(real, 0, Nyx::AxRe, 1, real.nGrow(), Ax_new.nGrow(), geom.periodicity(),FabArrayBase::COPY);
    Ax_new.ParallelCopy(imag, 0, Nyx::AxIm, 1, imag.nGrow(), Ax_new.nGrow(), geom.periodicity(),FabArrayBase::COPY);

    for (MFIter mfi(Ax_new); mfi.isValid(); ++mfi)
      {
	FArrayBox& fab = Ax_new[mfi];
	Real* ax = fab.dataPtr();
	const Box& axbox = fab.box();
	fort_ax_fields(ax, axbox.loVect(), axbox.hiVect());
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

void Nyx::init_rhs(MultiFab& rhs, Geometry& geom)
{
    const Real* dx = geom.CellSize();
    Box domain(geom.Domain());
    int i = 0 ;

    int nc = rhs.nComp();

    for (MFIter mfi(rhs,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        fort_init_rhs(BL_TO_FORTRAN_BOX(tbx),
                      BL_TO_FORTRAN_ANYD(rhs[mfi]),
                      BL_TO_FORTRAN_BOX(domain),
                      geom.ProbLo(), geom.ProbHi(), dx);
        i++;
    }

    Real sum_rhs = rhs.sum();
    amrex::Print() << "Sum of rhs over the domain was " << sum_rhs << std::endl;

    sum_rhs = sum_rhs / domain.numPts();
    rhs.plus(-sum_rhs,0,1,0);

    sum_rhs = rhs.sum();
    amrex::Print() << "Sum of rhs over the domain is now " << sum_rhs << std::endl;
}

#endif //FDM
