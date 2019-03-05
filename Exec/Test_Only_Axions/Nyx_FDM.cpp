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
    //testing swfft
    amrex::IntVect n_cell;
    amrex::IntVect max_grid_size;
    amrex::MultiFab rhs;
    amrex::MultiFab lhs;
    amrex::Geometry geom;
    int verbose_sw=2;


    {
        ParmParse pp;

        // Read in n_cell.  Use defaults if not explicitly defined.
        int cnt = pp.countval("amr.n_cell");

        if (cnt > 1) {
            Vector<int> ncs;
            pp.getarr("amr.n_cell",ncs);
            n_cell = IntVect{ncs[0],ncs[1],ncs[2]};
        } else if (cnt > 0) {
            int ncs;
            pp.get("amr.n_cell",ncs);
            n_cell = IntVect{ncs,ncs,ncs};
        } else {
            std::cout << "WARNING: amr.n_cell not found in inputs, setting n_cell to 32^3! " << std::endl;
            n_cell = IntVect{32,32,32};
        }

        // Read in max_grid_size.  Use defaults if not explicitly defined.
        cnt = pp.countval("amr.max_grid_size");
        if (cnt > 1) {
            Vector<int> mgs;
            pp.getarr("amr.max_grid_size",mgs);
            max_grid_size = IntVect{mgs[0],mgs[1],mgs[2]};
        } else if (cnt > 0) {
            int mgs;
            pp.get("amr.max_grid_size",mgs);
            max_grid_size = IntVect{mgs,mgs,mgs};
        } else {
            std::cout << "WARNING: amr.max_grid_size not found in inputs, setting max_grid_size to 32^3! " << std::endl;
            max_grid_size = IntVect{32,32,32};
        }

        pp.query("verbose", verbose);
    }

    BoxArray ba;

    Real dx = 1./double(n_cell[0]);
    IntVect dom_lo(0,0,0);
    IntVect dom_hi(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1);
    Box domain(dom_lo,dom_hi);
    ba.define(domain);
    ba.maxSize(max_grid_size);
    Real x_hi = n_cell[0]*dx;
    Real y_hi = n_cell[1]*dx;
    Real z_hi = n_cell[2]*dx;
    RealBox real_box({0.0,0.0,0.0}, {x_hi,y_hi,z_hi});
    // The FFT assumes fully periodic boundaries
    std::array<int,3> is_periodic {1,1,1};
    geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());

    // Make sure we define both the soln and the rhs with the same DistributionMapping
    DistributionMapping dmap{ba};

    rhs.define(ba, dmap, 1, 0);
    lhs.define(ba, dmap, 1, 0);
    init_rhs(rhs, geom);

    // //Uncomment if you want to plot rhs before swfft
    // const std::string pltfile0 = "MFAB_plot000";
    // std::string componentName = "rhs";
    // writeMultiFabAsPlotFile(pltfile0, rhs, componentName);

    swfft_test(rhs, lhs, geom, verbose_sw);

    // //Uncomment if you want to plot lhs after swfft
    // const std::string pltfile1 = "MFAB_plot002";
    // std::string componentName1 = "rhs";
    // writeMultiFabAsPlotFile(pltfile1, lhs, componentName1);
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
