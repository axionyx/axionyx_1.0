#ifdef FDM
 
#include "Nyx.H"
#include "Nyx_F.H"
//#include <AMReX_Particles_F.H>
#include <AMReX_MultiFab.H>
#ifdef GRAVITY
#include "Gravity.H"
#include <Gravity_F.H>
#endif

void
Nyx::advance_FDM_FD (amrex::Real time,
                      amrex::Real dt,
                      amrex::Real a_old,
                      amrex::Real a_new)
{
    // amrex::Real se, ske;
    // const amrex::Real prev_time    = state[State_Type].prevTime();
    // const amrex::Real cur_time     = state[State_Type].curTime();
    // const int  finest_level = parent->finestLevel();

    const amrex::Real a_half = 0.5 * (a_old + a_new);

    amrex::MultiFab&  Ax_old = get_old_data(Axion_Type);
    amrex::MultiFab&  Ax_new = get_new_data(Axion_Type);
    amrex::MultiFab& Phi_old = get_old_data(PhiGrav_Type);
    amrex::MultiFab& Phi_new = get_new_data(PhiGrav_Type);

    if (verbose && amrex::ParallelDescriptor::IOProcessor() ){
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
    if (Phi_old.contains_nan(0, 1, 0))
      amrex::Abort("Phi_new has NaNs::advance_FDM_FD()");
#endif

    const amrex::Real* dx      = geom.CellSize();
    amrex::Real        courno  = -1.0e+200;
    //Need the positions of the physical boundaries for the 'sponge' neccessary for apropriate 
    //boundary conditions for the Schroedinger-Poisson system (arXiv:gr-qc/0404014v2)
    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* prob_hi = geom.ProbHi();


    // Define the gravity vector so we can pass this to ca_umdrv.
    // MultiFab grav_vector(grids, BL_SPACEDIM, 3, Fab_allocate);
    const auto& dm = get_level(level).get_new_data(Axion_Type).DistributionMap();
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
	     fpi(*this,  Ax_new, 4, time, Axion_Type,   0, Nyx::NUM_AX),
	     pfpi(*this, Phi_old, 4, time, PhiGrav_Type, 0, 1);
	     fpi.isValid() && pfpi.isValid();
	     ++fpi,++pfpi)
       {
        // Create FAB for extended grid values (including boundaries) and fill.
        // const int  mfi_index = fpi.index();
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
             // BL_TO_FORTRAN(grav_vector[fpi]),
             BL_TO_FORTRAN(phiold),
             dx, prob_lo, prob_hi, &dt,
             &cflLoc, &a_old, &a_half, &a_new, verbose);
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

#endif //FDM
