#ifdef FDM
 
#include "Nyx.H"
#include "Nyx_F.H"
#include <AMReX_Particles_F.H>
using namespace amrex;


#ifdef GRAVITY
#	include "Gravity.H"
#endif
 
using std::string;

Real
Nyx::advance_FDM (Real time,
		  Real dt,
		  int  iteration,
		  int  ncycle)
{
  BL_PROFILE("Nyx::advance_particles_only()");

  // A particle in cell (i) can affect cell values in (i-1) to (i+1)                                                                                                                                             
  int stencil_deposition_width = 1;

  // For FDM Gaussian kernels this is increased to
  int stencil_deposition_width_fdm = ceil(sigma_ax*theta_ax)*pow(2,level);

  // A particle in cell (i) may need information from cell values in (i-1) to (i+1)                                                                                                                              
  //   to update its position (typically via interpolation of the acceleration from the grid)                                                                                                                    
  int stencil_interpolation_width = 1;

  // A particle that starts in cell (i + ncycle) can reach                                                                                                                                                       
  //   cell (i) in ncycle number of steps .. after "iteration" steps                                                                                                                                             
  //   the particle has to be within (i + ncycle+1-iteration) to reach cell (i)                                                                                                                                  
  //   in the remaining (ncycle-iteration) steps                                                                                                                                                                 

  // *** ghost_width ***  is used                                                                                                                                                                                
  //   *) to set how many cells are used to hold ghost particles i.e copies of particles                                                                                                                         
  //      that live on (level-1) can affect the grid over all of the ncycle steps.                                                                                                                               
  //      We define ghost cells at the coarser level to cover all iterations so                                                                                                                                  
  //      we can't reduce this number as iteration increases.                                                                                                                                                    

  // ghost_width for FDM Gaussian kernels is separately set in setup_ghost_particles

  int ghost_width = ncycle + stencil_deposition_width;

  // *** where_width ***  is used                                                                                                                                                                                
  //   *) to set how many cells the Where call in moveKickDrift tests =                                                                                                                                          
  //      ghost_width + (1-iteration) - 1:                                                                                                                                                                       
  //      the minus 1 arises because this occurs *after* the move                                                                                                                                                

  int where_width =  ghost_width + (1-iteration)  - 1;

  // *** grav_n_grow *** is used                                                                                                                                                                                 
  //   *) to determine how many ghost cells we need to fill in the MultiFab from                                                                                                                                 
  //      which the particle interpolates its acceleration                                                                                                                                                       
  //   *) to set how many cells the Where call in moveKickDrift tests = (grav.nGrow()-2).                                                                                                                        

  int grav_n_grow = ncycle+stencil_deposition_width_fdm + (1-iteration) +
    stencil_interpolation_width ;

    bool show_timings=false;
    //I put this here; its value is set differently in Nyx_advance (advance_hydro_plus_particles), but in the old version of the code, it was set to 3.
    // grav_n_grow = 3;
    // Sanity checks
    if (do_hydro)
        amrex::Abort("In `advance_particles_only` but `do_hydro` is true");

    const int finest_level = parent->finestLevel();
    int finest_level_to_advance;
    bool nosub = !parent->subCycle();
    amrex::Real dt_lev;
    const amrex::Real strt = ParallelDescriptor::second();

    if (nosub)
    {
        if (level > 0)
            return dt;
            
        finest_level_to_advance = finest_level;
    }
    else
    {
        if (strict_subcycling)
        {
            finest_level_to_advance = level;
        }
        else
        {
            // This level was advanced by a previous multilevel advance.
            if (level > 0 && ncycle == 1)
                return dt;
            
            // Find the finest level to advance
            int lev = level;
            while(lev < finest_level && parent->nCycle(lev+1) == 1)
                lev++;
            finest_level_to_advance = lev;
        }
//TODO; check if we need the stuff below. I commented it, as it concerns dealing with nbody particles. The axionic part is below.
        // We must setup virtual and Ghost Particles
        //
        // Setup the virtual particles that represent finer level particles
        // 
        setup_virtual_particles();
        //
        // Setup ghost particles for use in finer levels. Note that Ghost particles
        // that will be used by this level have already been created, the
        // particles being set here are only used by finer levels.
        //
       for(int lev = level; lev <= finest_level_to_advance && lev < finest_level; lev++)
       {
          get_level(lev).setup_ghost_particles(ghost_width);
       }
    }



    //
    // Move current data to previous, clear current.
    // Don't do this if a coarser level has done this already.
    //
    if (level == 0 || iteration > 1)
    {
        for (int lev = level; lev <= finest_level; lev++)
        {
            dt_lev = parent->dtLevel(lev);
            for (int k = 0; k < NUM_STATE_TYPE; k++)
            {
                get_level(lev).state[k].allocOldData();
                get_level(lev).state[k].swapTimeLevels(dt_lev);
            }

	    //If we don't evolve hydro states we simple copy them
	    if(!do_hydro){
	      MultiFab& S_old = get_level(lev).get_old_data(State_Type);
	      MultiFab& S_new = get_level(lev).get_new_data(State_Type);
	      MultiFab::Copy(S_new, S_old, 0, 0, S_old.nComp(), 0);
	      
	      MultiFab& D_old = get_level(lev).get_old_data(DiagEOS_Type);
	      MultiFab& D_new = get_level(lev).get_new_data(DiagEOS_Type);
	      MultiFab::Copy(D_new, D_old, 0, 0, D_old.nComp(), 0);
	    }
        }
    }

    const amrex::Real prev_time = state[State_Type].prevTime();
    const amrex::Real cur_time  = state[State_Type].curTime();

    const amrex::Real a_old     = get_comoving_a(prev_time);
    const amrex::Real a_new     = get_comoving_a(cur_time);

//TODO: once again, this part only concerns particles.
#ifdef GRAVITY
    //
    // We now do a multilevel solve for old Gravity. This goes to the 
    // finest level regardless of subcycling behavior. Consequentially,
    // If we are subcycling we skip this step on the first iteration of
    // finer levels.
    if (level == 0 || iteration > 1)
    {
        // fix fluxes on finer grids
        if (do_reflux)
        {
            for (int lev = level; lev < finest_level; lev++)
            {
                gravity->zero_phi_flux_reg(lev + 1);
            }
        }
        
        // swap grav data
        for (int lev = level; lev <= finest_level; lev++)
            get_level(lev).gravity->swap_time_levels(lev);

        if (show_timings)
        {
            const int IOProc = ParallelDescriptor::IOProcessorNumber();
            amrex::Real end = ParallelDescriptor::second() - strt;
            ParallelDescriptor::ReduceRealMax(end,IOProc);
            if (ParallelDescriptor::IOProcessor())
               std::cout << "Time before solve for old phi " << end << '\n';
        }

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "\n... old-time level solve at level " << level << '\n';
            
        //
        // Solve for phi
        // If a single-level calculation we can still use the previous phi as a guess.
        // TODO: Check this.
        int use_previous_phi_as_guess = 0;
        gravity->multilevel_solve_for_old_phi(level, finest_level,
                                              use_previous_phi_as_guess);
    }
   //
   // Advance Particles
   //
   if (Nyx::theActiveParticles().size() > 0)
   {
       // Advance the particle velocities to the half-time and the positions to the new time
       // We use the cell-centered gravity to correctly interpolate onto particle locations
       if (particle_move_type == "Gravitational")
       {
           const amrex::Real a_half = 0.5 * (a_old + a_new);

           if (particle_verbose && ParallelDescriptor::IOProcessor())
               std::cout << "moveKickDrift ... updating particle positions and velocity\n";

           for (int lev = level; lev <= finest_level_to_advance; lev++)
           {
               // We need grav_n_grow grow cells to track boundary particles
               const auto& ba = get_level(lev).get_new_data(State_Type).boxArray();
               const auto& dm = get_level(lev).get_new_data(PhiGrav_Type).DistributionMap();
               MultiFab grav_vec_old(ba, dm, BL_SPACEDIM, grav_n_grow);
               get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old, time);
	       MultiFab& phi = get_level(lev).get_old_data(PhiGrav_Type);

               for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
                   Nyx::theActiveParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half);
	       Nyx::theFDMPC()->moveKickDriftFDM(phi, grav_n_grow, grav_vec_old, lev, dt, a_old, a_half);

               // // Only need the coarsest virtual particles here.
               // if (lev == level && level < finest_level){
                   for (int i = 0; i < Nyx::theVirtualParticles().size(); i++)
                       Nyx::theVirtualParticles()[i]->moveKickDrift(grav_vec_old, level, dt, a_old, a_half);
		   Nyx::theGhostFDMPC()->moveKickDriftFDM(phi, grav_n_grow, grav_vec_old, lev, dt, a_old, a_half);
	       // }

               // Miiiight need all Ghosts
               for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
                   Nyx::theGhostParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half);
	       Nyx::theGhostFDMPC()->moveKickDriftFDM(phi, grav_n_grow, grav_vec_old, lev, dt, a_old, a_half);
           }
       }
   }

#endif

#ifndef FDM_GB
    //Advance Axions
    for (int lev = level; lev <= finest_level_to_advance; lev++)
        get_level(lev).advance_FDM_FD(time, dt, a_old, a_new);

    // Always average down from finer to coarser.
    for (int lev = finest_level_to_advance-1; lev >= level; lev--)
        get_level(lev).average_down(Axion_Type);
#endif

    //                                                                                                                                                                                                           
    // Call the hydro advance at each level to be advanced                                                                                                                                                       
    //                                                                                                                                                                                                            
    if(do_hydro){
      BL_PROFILE_VAR("just_the_hydro", just_the_hydro);
      for (int lev = level; lev <= finest_level_to_advance; lev++)
	{
#ifdef SDC
	  if (sdc_split > 0)
	    {
	      get_level(lev).sdc_hydro(time, dt, a_old, a_new);
	    } else {
	    get_level(lev).strang_hydro(time, dt, a_old, a_new);
	  }
#else
	  get_level(lev).strang_hydro(time, dt, a_old, a_new);
#endif
	}
      BL_PROFILE_VAR_STOP(just_the_hydro);
      
      //                                                                                                                                                                                                        
      // We must reflux before doing the next gravity solve                                                                                                                                                     
      //                                                                                                                                                                                                        
      if (do_reflux)
	{
	  for (int lev = level; lev < finest_level_to_advance; lev++)
	    {
	      get_level(lev).reflux();
	    }
	}
      
      // Always average down the new state from finer to coarser.                                                                                                                                              
      for (int lev = finest_level_to_advance-1; lev >= level; lev--)
	{
	  get_level(lev).average_down(  State_Type);
	  get_level(lev).average_down(DiagEOS_Type);
	}
    }



#ifdef GRAVITY

    //
    // Here we use the "old" phi from the current time step as a guess for this
    // solve
    //
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
        MultiFab::Copy(parent->getLevel(lev).get_new_data(PhiGrav_Type),
                       parent->getLevel(lev).get_old_data(PhiGrav_Type),
                       0, 0, 1, parent->getLevel(lev).get_old_data(PhiGrav_Type).nGrow());
    }

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Time before solve for new phi " << end << '\n';
    }

    // Solve for new Gravity                                                                                                                                                                                    
    BL_PROFILE_VAR("solve_for_new_phi", solve_for_new_phi);
    int use_previous_phi_as_guess = 1;
    if (finest_level_to_advance > level)
      {
        // The particle may be as many as "iteration" ghost cells out                                                                                                                                         
        int ngrow_for_solve = iteration + stencil_deposition_width;
        gravity->multilevel_solve_for_new_phi(level, finest_level_to_advance,
                                              ngrow_for_solve,
                                              use_previous_phi_as_guess);
      }
    else
      {
        int fill_interior = 0;
        gravity->solve_for_new_phi(level,get_new_data(PhiGrav_Type),
				   gravity->get_grad_phi_curr(level),
				   fill_interior, grav_n_grow);
      }
    BL_PROFILE_VAR_STOP(solve_for_new_phi);


    if(do_hydro){
      // Reflux                                                                                                                                                                                                 
      if (do_reflux)
	{
	  for (int lev = level; lev <= finest_level_to_advance; lev++)
	    {
	      gravity->add_to_fluxes(lev, iteration, ncycle);
	    }
	}
      
      //                                                                                                                                                                                                       
      // Now do corrector part of source term update                                                                                                                                                           
      //                                                                                                                                                                                                       
      for (int lev = level; lev <= finest_level_to_advance; lev++)
	{
	  MultiFab& S_old = get_level(lev).get_old_data(State_Type);
	  MultiFab& S_new = get_level(lev).get_new_data(State_Type);
	  MultiFab& D_new = get_level(lev).get_new_data(DiagEOS_Type);
	  MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
	  reset_e_src.setVal(0.0);
	  
	  const auto& ba = get_level(lev).get_new_data(State_Type).boxArray();
	  const auto& dm = get_level(lev).get_new_data(State_Type).DistributionMap();
	  
	  // These vectors are only used for the call to correct_gsrc so they                                                                                                                                     
	  //    don't need any ghost cells                                                                                                                                                                      
	  MultiFab grav_vec_old(ba, dm, BL_SPACEDIM, 0);
	  MultiFab grav_vec_new(ba, dm, BL_SPACEDIM, 0);
	  
	  get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old, time);
	  get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new, cur_time);
	  
	  const Real* dx = get_level(lev).Geom().CellSize();
	  
#ifdef _OPENMP
#pragma omp parallel
#endif
	  for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	    {
	      const Box& bx = mfi.tilebox();
	      
	      Real se  = 0;
	      Real ske = 0;
	      
	      fort_correct_gsrc
		(bx.loVect(), bx.hiVect(), BL_TO_FORTRAN(grav_vec_old[mfi]),
		 BL_TO_FORTRAN(grav_vec_new[mfi]), BL_TO_FORTRAN(S_old[mfi]),
		 BL_TO_FORTRAN(S_new[mfi]), &a_old, &a_new, &dt);
	    }
	  
	  // First reset internal energy before call to compute_temp                                                                                                                                        
	  get_level(lev).reset_internal_energy(S_new,D_new,reset_e_src);
	  get_level(lev).compute_new_temp(S_new,D_new);
	}

      // TODO: same here
      
      // Must average down again after doing the gravity correction;
      //      always average down from finer to coarser.
      for (int lev = finest_level_to_advance-1; lev >= level; lev--)
	get_level(lev).average_down();  
    }

    if (Nyx::theActiveParticles().size() > 0)
      {
	// Advance the particle velocities by dt/2 to the new time. We use the
	// cell-centered gravity to correctly interpolate onto particle
	// locations.
	if (particle_move_type == "Gravitational")
	  {
	    const amrex::Real a_half = 0.5 * (a_old + a_new);
	    
	    if (particle_verbose && ParallelDescriptor::IOProcessor())
	      std::cout << "moveKick ... updating velocity only\n";
	    
	    for (int lev = level; lev <= finest_level_to_advance; lev++)
	      {
		const BoxArray& ba = get_level(lev).get_new_data(State_Type).boxArray();
		const auto& dm = get_level(lev).get_new_data(PhiGrav_Type).DistributionMap();
		MultiFab grav_vec_new(ba, dm, BL_SPACEDIM, grav_n_grow);
		get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new, cur_time);
		MultiFab& phi = get_level(lev).get_new_data(PhiGrav_Type);
		
		for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
		  Nyx::theActiveParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);
		Nyx::theFDMPC()->moveKickFDM(phi, grav_n_grow, grav_vec_new, lev, dt, a_new, a_half);
		
		// // Only need the coarsest virtual particles here.
		// if (lev == level && level < finest_level){
		  for (int i = 0; i < Nyx::theVirtualParticles().size(); i++)
		    Nyx::theVirtualParticles()[i]->moveKick(grav_vec_new, level, dt, a_new, a_half);
		  Nyx::theGhostFDMPC()->moveKickFDM(phi, grav_n_grow, grav_vec_new, lev, dt, a_new, a_half);
		// }
		
		// Miiiight need all Ghosts
		for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
		  Nyx::theGhostParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);
		Nyx::theGhostFDMPC()->moveKickFDM(phi, grav_n_grow, grav_vec_new, lev, dt, a_new, a_half);

	      }
	  }
      }
    
    if(!do_hydro){
      //we need to get new grids, since it implicitely fills the new part of the
      //Gravity_Type, which we need for the interpolated values of G at higher 
      //levels when filling ghosts from lower levels.
      const auto& dm = get_level(level).get_new_data(State_Type).DistributionMap();
      MultiFab grav_vec_new(grids, dm, BL_SPACEDIM, grav_n_grow);
      get_level(level).gravity->get_new_grav_vector(level, grav_vec_new, cur_time);
    }

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Time  after moveKick " << end << '\n';
    }
#endif

#ifdef FDM_GB
    for (int lev = level; lev <= finest_level_to_advance; lev++){
    //Define neccessary number of ghost cells                                                                                                                                                                   
    int ng = ceil(Nyx::sigma_ax*Nyx::theta_ax)*pow(2,lev);

    //Initialize MultiFabs                                                                                                                                                                                     
    MultiFab& Ax_new = get_level(lev).get_new_data(Axion_Type);
    Ax_new.setVal(0.);
    MultiFab fdmreal(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
    fdmreal.setVal(0.);
    MultiFab fdmimag(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
    fdmimag.setVal(0.);

    //Deposit Gaussian Beams                                                                                                                                                                               
    if(Nyx::theFDMPC())
      Nyx::theFDMPC()->DepositFDMParticles(fdmreal,fdmimag,lev);
    if(Nyx::theGhostFDMPC())
      Nyx::theGhostFDMPC()->DepositFDMParticles(fdmreal,fdmimag,lev);
    if(Nyx::theVirtFDMPC())
      Nyx::theVirtFDMPC()->DepositFDMParticles(fdmreal,fdmimag,lev);

    //Update real part in FDM state                                                                                                                                                                             
    Ax_new.ParallelCopy(fdmreal, 0, Nyx::AxRe, 1, fdmreal.nGrow(),
                        Ax_new.nGrow(), parent->Geom(lev).periodicity(),FabArrayBase::ADD);

    //Update imaginary part in FDM state                                                                                                                                                                        
    Ax_new.ParallelCopy(fdmimag, 0, Nyx::AxIm, 1, fdmimag.nGrow(),
                        Ax_new.nGrow(), parent->Geom(lev).periodicity(),FabArrayBase::ADD);

    //Update density in FDM state                                                                                                                                                                              
    AmrLevel* amrlev = &parent->getLevel(lev);
    for (amrex::FillPatchIterator fpi(*amrlev,  Ax_new); fpi.isValid(); ++fpi)
      {
        if (Ax_new[fpi].contains_nan())
	  amrex::Abort("Nans in state just before FDM density update");
        BL_FORT_PROC_CALL(FORT_FDM_FIELDS, fort_fdm_fields)
          (BL_TO_FORTRAN(Ax_new[fpi]));
        if (Ax_new[fpi].contains_nan())
	  amrex::Abort("Nans in state just after FDM density update");
      }
    }
#endif

    if(do_hydro){
      //                                                                                                                                                                                                       
      // Synchronize Energies                                                                                                                                                                                  
      //                                                                                                                                                                                                        
      for (int lev = level; lev <= finest_level_to_advance; lev++)
	{
	  MultiFab& S_new = get_level(lev).get_new_data(State_Type);
	  MultiFab& D_new = get_level(lev).get_new_data(DiagEOS_Type);
	  MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
	  reset_e_src.setVal(0.0);
	  
	  get_level(lev).reset_internal_energy(S_new,D_new,reset_e_src);
	  
	}
    }

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Time  at end of routine " << end << '\n';
    }

    // Redistribution happens in post_timestep
    return dt;
}
#endif
