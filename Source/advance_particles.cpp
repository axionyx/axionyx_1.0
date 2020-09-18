#ifdef GRAVITY
 
#include <AMReX_BLProfiler.H>
#include "Nyx.H"
#include "Nyx_F.H"
#include "Gravity.H"
#include <Gravity_F.H>

using namespace amrex;
 
using std::string;

Real
Nyx::advance_particles_only (Real time,
                             Real dt,
                             int  iteration,
                             int  ncycle)

  // Arguments:
  //    time      : the current simulation time
  //    dt        : the timestep to advance (e.g., go from time to
  //                time + dt)
  //    iteration : where we are in the current AMR subcycle.  Each
  //                level will take a number of steps to reach the
  //                final time of the coarser level below it.  This
  //                counter starts at 1
  //    ncycle    : the number of subcycles at this level

{
     BL_PROFILE("Nyx::advance_particles_only()");

    // A particle in cell (i) can affect cell values in (i-1) to (i+1)
    // ! Now in Nyx.H
    // int stencil_deposition_width = 1;

// #ifdef FDM
//     // For FDM Gaussian kernels this is increased to                                                                                                                                                          
//     int stencil_deposition_width_fdm = ceil(sigma_fdm*theta_fdm/get_level(level).Geom().CellSize()[0]);
// #endif 

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
 
    int ghost_width = ncycle + stencil_deposition_width;
    if(level<parent->finestLevel())
      ghost_width = parent->nCycle(level+1) + stencil_deposition_width;

    // *** where_width ***  is used
    //   *) to set how many cells the Where call in moveKickDrift tests =
    //      ghost_width + (1-iteration) - 1:
    //      the minus 1 arises because this occurs *after* the move

    int where_width =  ghost_width + (1-iteration)  - 1;
 
    // *** grav_n_grow *** is used
    //   *) to determine how many ghost cells we need to fill in the MultiFab from
    //      which the particle interpolates its acceleration
    //   *) to set how many cells the Where call in moveKickDrift tests = (grav.nGrow()-2).
    //   *) the (1-iteration) arises because the ghost particles are created on the coarser
    //      level which means in iteration 2 the ghost particles may have moved 1 additional cell along
 
    int grav_n_grow = ghost_width + (1-iteration) + (iteration-1) +
                      stencil_interpolation_width ;
#ifdef FDM
    // Plus one since we need to take the derivative of grav_vector                                                                                                                                            
    // in order to obtain hession of potential
    grav_n_grow += 1;
#endif

    // Sanity checks
    if (do_hydro)
        amrex::Abort("In `advance_particles_only` but `do_hydro` is true");

    if (!do_grav)
        amrex::Abort("In `advance_particles_only` but `do_grav` not true");
    const int finest_level = parent->finestLevel();
    int finest_level_to_advance;
    bool nosub = !parent->subCycle();
    
    if (nosub)
    {
        if (level > 0)
            return dt;
            
        finest_level_to_advance = finest_level;

#ifdef FDM
	// Need virtual and ghost particles for Gaussian beam deposition
        for(int lev = level; lev < finest_level; lev++){
	  if(levelmethod[lev]==GBlevel || levelmethod[lev]==CWlevel)
	    setup_virtual_particles();
	  if(levelmethod[lev+1]==GBlevel || levelmethod[lev]==CWlevel)
	    get_level(lev).setup_ghost_particles(ghost_width);
	}
#endif

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

    Real dt_lev;

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
#ifndef NO_HYDRO
            MultiFab& S_old = get_level(lev).get_old_data(State_Type);
            MultiFab& S_new = get_level(lev).get_new_data(State_Type);
            MultiFab::Copy(S_new, S_old, 0, 0, S_old.nComp(), 0);
#endif
        }
    }

    const Real prev_time = state[Gravity_Type].prevTime();
    const Real cur_time  = state[Gravity_Type].curTime();

    const Real a_old     = get_comoving_a(prev_time);
    const Real a_new     = get_comoving_a(cur_time);

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

        //
        // Solve for phi
        // If a single-level calculation we can still use the previous phi as a guess.
        // TODO: Check this.
        // int use_previous_phi_as_guess = 1;
        int use_previous_phi_as_guess = 0;

#ifdef CGRAV
        MultiFab& phi_old = get_level(level).get_old_data(PhiGrav_Type);
        MultiFab& phi_new = get_level(level).get_new_data(PhiGrav_Type);

        prescribe_grav_potential(phi_old, Geom() , level, finest_level);

	MultiFab::Copy(phi_new, phi_old, 0, 0, phi_old.nComp(), 0);


#else

// #ifdef FDM
//         // for(int lev = level; lev <= finest_level_to_advance && lev < finest_level; lev++)
//         for(int lev = finest_level_to_advance+1; lev < finest_level; lev++)
//         {
// 	  if(levelmethod[lev+1]==Nyx::GBlevel)
// 	    get_level(lev).setup_ghost_particles(ghost_width);
//         }
// #endif

        gravity->multilevel_solve_for_old_phi(level, finest_level, grav_n_grow,
                                              use_previous_phi_as_guess);

// #ifdef FDM
//         // for(int lev = level; lev <= finest_level_to_advance && lev < finest_level; lev++)
//         for(int lev = finest_level_to_advance+1; lev < finest_level; lev++)
//         {
// 	  if(levelmethod[lev+1]==Nyx::GBlevel){
// 	    if(theGhostFDMPC())
// 	      theGhostFDMPC()->RemoveParticlesAtLevel(lev+1);
// 	  if(theGhostFDMwkbPC())
// 	    theGhostFDMwkbPC()->RemoveParticlesAtLevel(lev+1);
// 	  }
// 	}
// #endif

#endif
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
            const Real a_half = 0.5 * (a_old + a_new);

            if (particle_verbose && ParallelDescriptor::IOProcessor())
                std::cout << "moveKickDrift ... updating particle positions and velocity\n";

            for (int lev = level; lev <= finest_level_to_advance; lev++)
            {
                // We need grav_n_grow grow cells to track boundary particles
                const auto& ba = get_level(lev).get_old_data(Gravity_Type).boxArray();
                const auto& dm = get_level(lev).get_old_data(Gravity_Type).DistributionMap();
                MultiFab grav_vec_old(ba, dm, BL_SPACEDIM, grav_n_grow);
                get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old, prev_time);
                
                for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
                    Nyx::theActiveParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half, iteration);

                // Only need the coarsest virtual particles here.
                if (lev == level && level < finest_level)
                    for (int i = 0; i < Nyx::theVirtualParticles().size(); i++)
                        Nyx::theVirtualParticles()[i]->moveKickDrift(grav_vec_old, level, dt, a_old, a_half,iteration);

                // Miiiight need all Ghosts
                for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
                    Nyx::theGhostParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_new, a_half,where_width);

#ifdef FDM
		MultiFab& Phi_old = get_level(lev).get_old_data(PhiGrav_Type);
		if(Nyx::theFDMPC())
		  Nyx::theFDMPC()->moveKickDriftFDM(Phi_old, grav_n_grow, grav_vec_old, lev, dt, a_old, a_half,iteration);
		if(Nyx::theFDMwkbPC())
		  Nyx::theFDMwkbPC()->moveKickDriftFDM(Phi_old, grav_n_grow, grav_vec_old, lev, dt, a_old, a_half,iteration);
		if(Nyx::theFDMphasePC())
		  Nyx::theFDMphasePC()->moveKickDriftFDM(Phi_old, grav_n_grow, grav_vec_old, lev, dt, a_old, a_half,iteration);

		//Need to do this only when reconstructing the density from GBs on this level.
		if(levelmethod[lev]==GBlevel){
		  int ghost_width_fdm = parent->nCycle(lev)+ceil(Nyx::sigma_fdm*Nyx::theta_fdm/get_level(lev).Geom().CellSize()[0]);
		  int where_width_fdm =  ghost_width_fdm + (1-iteration)  - 1;
		  int grav_n_grow_fdm = ghost_width_fdm + stencil_interpolation_width + 1;
		  MultiFab grav_vec_old_fdm(ba, dm, BL_SPACEDIM, grav_n_grow_fdm);
		  get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old_fdm, prev_time);

		  if(Nyx::theGhostFDMPC())
		    Nyx::theGhostFDMPC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,where_width_fdm);
		  if(Nyx::theVirtFDMPC())
		    Nyx::theVirtFDMPC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,iteration);
		  if(Nyx::theGhostFDMwkbPC())
		    Nyx::theGhostFDMwkbPC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,where_width_fdm);
		  if(Nyx::theVirtFDMwkbPC())
		    Nyx::theVirtFDMwkbPC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,iteration);
		  if(Nyx::theGhostFDMphasePC())
		    Nyx::theGhostFDMphasePC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,where_width_fdm);
		  if(Nyx::theVirtFDMphasePC())
		    Nyx::theVirtFDMphasePC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,iteration);
		}
		else if (levelmethod[lev]==CWlevel){
		  int ghost_width_fdm = parent->nCycle(lev)+ceil(Nyx::theta_fdm);
		  int where_width_fdm =  ghost_width_fdm + (1-iteration)  - 1;
		  int grav_n_grow_fdm = ghost_width_fdm + stencil_interpolation_width + 1;
		  MultiFab grav_vec_old_fdm(ba, dm, BL_SPACEDIM, grav_n_grow_fdm);
		  get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old_fdm, prev_time);

		  if(Nyx::theGhostFDMPC())
		    Nyx::theGhostFDMPC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,where_width_fdm);
		  if(Nyx::theVirtFDMPC())
		    Nyx::theVirtFDMPC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,iteration);
		  if(Nyx::theGhostFDMwkbPC())
		    Nyx::theGhostFDMwkbPC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,where_width_fdm);
		  if(Nyx::theVirtFDMwkbPC())
		    Nyx::theVirtFDMwkbPC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,iteration);
		  if(Nyx::theGhostFDMphasePC())
		    Nyx::theGhostFDMphasePC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,where_width_fdm);
		  if(Nyx::theVirtFDMphasePC())
		    Nyx::theVirtFDMphasePC()->moveKickDriftFDM(Phi_old, grav_n_grow_fdm, grav_vec_old_fdm, lev, dt, a_old, a_half,iteration);
		}
#endif
            }
        }
    }

#ifdef FDM
    //Advance FDM wavefunction with finite difference method on FDlevels                                                                                                        
    for (int lev = level; lev <= finest_level_to_advance; lev++)
      if(levelmethod[lev]==FDlevel)
	get_level(lev).advance_FDM_FD(time, dt, a_old, a_new);
      else if(levelmethod[lev]==PSlevel)
	get_level(lev).advance_FDM_PS(time, dt, a_old, a_new);

    // Always average down from finer to coarser.                                                                                                                                                
    for (int lev = finest_level_to_advance-1; lev >= level; lev--)
      get_level(lev).average_down(Axion_Type);
#endif

    //
    // Here we use the "old" phi from the current time step as a guess for this
    // solve
    //
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
        MultiFab::Copy(parent->getLevel(lev).get_new_data(PhiGrav_Type),
                       parent->getLevel(lev).get_old_data(PhiGrav_Type),
                       0, 0, 1, 0);
    }

    // Solve for new Gravity
    int use_previous_phi_as_guess = 1;
    if (finest_level_to_advance > level)
    {
      gravity->multilevel_solve_for_new_phi(level, finest_level_to_advance, grav_n_grow,
                                              use_previous_phi_as_guess);
    }
    else
    {
        int fill_interior = 0;
        gravity->solve_for_new_phi(level,get_new_data(PhiGrav_Type),
                               gravity->get_grad_phi_curr(level),
                               fill_interior, grav_n_grow);
    }

    if (Nyx::theActiveParticles().size() > 0)
    {
        // Advance the particle velocities by dt/2 to the new time. We use the
        // cell-centered gravity to correctly interpolate onto particle
        // locations.
        if (particle_move_type == "Gravitational")
        {
            const Real a_half = 0.5 * (a_old + a_new);

            if (particle_verbose && ParallelDescriptor::IOProcessor())
                std::cout << "moveKick ... updating velocity only\n";

            for (int lev = level; lev <= finest_level_to_advance; lev++)
            {
                const auto& ba = get_level(lev).get_new_data(PhiGrav_Type).boxArray();
                const auto& dm = get_level(lev).get_new_data(PhiGrav_Type).DistributionMap();
                MultiFab grav_vec_new(ba, dm, BL_SPACEDIM, grav_n_grow);
                get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new, cur_time);

                for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
                    Nyx::theActiveParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);

                // Virtual particles will be recreated, so we need not kick them.

                // Ghost particles need to be kicked except during the final iteration.
                if (iteration != ncycle)
                    for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
                        Nyx::theGhostParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);

#ifdef FDM
		MultiFab& Phi_new = get_level(lev).get_new_data(PhiGrav_Type);
		if(Nyx::theFDMPC())
		  Nyx::theFDMPC()->moveKickFDM(Phi_new, grav_n_grow, grav_vec_new, lev, dt, a_new, a_half);
		if(Nyx::theFDMwkbPC())
		  Nyx::theFDMwkbPC()->moveKickFDM(Phi_new, grav_n_grow, grav_vec_new, lev, dt, a_new, a_half);
		if(Nyx::theFDMphasePC())
		  Nyx::theFDMphasePC()->moveKickFDM(Phi_new, grav_n_grow, grav_vec_new, lev, dt, a_new, a_half);

		//Need to do this only when reconstructing the density from GBs on this level.
		if(levelmethod[lev]==GBlevel){
		  int ghost_width_fdm = parent->nCycle(lev)+ceil(Nyx::sigma_fdm*Nyx::theta_fdm/get_level(lev).Geom().CellSize()[0]);
		  int where_width_fdm =  ghost_width_fdm + (1-iteration)  - 1;
		  int grav_n_grow_fdm = ghost_width_fdm + stencil_interpolation_width + 1;
		  MultiFab grav_vec_new_fdm(ba, dm, BL_SPACEDIM, grav_n_grow_fdm);
		  get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new_fdm, cur_time);

		  if(Nyx::theGhostFDMPC())
		    Nyx::theGhostFDMPC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theVirtFDMPC())
		    Nyx::theVirtFDMPC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theGhostFDMwkbPC())
		    Nyx::theGhostFDMwkbPC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theVirtFDMwkbPC())
		    Nyx::theVirtFDMwkbPC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theGhostFDMphasePC())
		    Nyx::theGhostFDMphasePC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theVirtFDMphasePC())
		    Nyx::theVirtFDMphasePC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		}
		else if (levelmethod[lev]==CWlevel){
		  int ghost_width_fdm = parent->nCycle(lev)+ceil(Nyx::theta_fdm);
		  int where_width_fdm =  ghost_width_fdm + (1-iteration)  - 1;
		  int grav_n_grow_fdm = ghost_width_fdm + stencil_interpolation_width + 1;
		  MultiFab grav_vec_new_fdm(ba, dm, BL_SPACEDIM, grav_n_grow_fdm);
		  get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new_fdm, cur_time);

		  if(Nyx::theGhostFDMPC())
		    Nyx::theGhostFDMPC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theVirtFDMPC())
		    Nyx::theVirtFDMPC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theGhostFDMwkbPC())
		    Nyx::theGhostFDMwkbPC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theVirtFDMwkbPC())
		    Nyx::theVirtFDMwkbPC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theGhostFDMphasePC())
		    Nyx::theGhostFDMphasePC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		  if(Nyx::theVirtFDMphasePC())
		    Nyx::theVirtFDMphasePC()->moveKickFDM(Phi_new, grav_n_grow_fdm, grav_vec_new_fdm, lev, dt, a_new, a_half);
		}
#endif

            }
        }
    }

#ifdef FDM
    for (int lev = level; lev <= finest_level_to_advance; lev++){

      //Only construct wavefunction from Gaussian beams on GBlevels
      if(levelmethod[lev]==GBlevel){

      amrex::Print() << "levelSteps " << lev << " "<< parent->levelSteps(lev) << '\n';

      //Define neccessary number of ghost cells                                                                                                                                                                 
      int ng = parent->nCycle(lev)+2.0*ceil(Nyx::sigma_fdm*Nyx::theta_fdm/get_level(lev).Geom().CellSize()[0]);
      

      //Initialize MultiFabs                                                                                                                                                                                    
      MultiFab& Ax_new = get_level(lev).get_new_data(Axion_Type);
      Ax_new.setVal(0.);
      MultiFab fdmreal(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
      fdmreal.setVal(0.);
      MultiFab fdmimag(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
      fdmimag.setVal(0.);

      //Deposit Gaussian Beams                                                                                                                                                                                   
      if(Nyx::theFDMPC())
      	Nyx::theFDMPC()->DepositFDMParticles(fdmreal,fdmimag,lev,a_new);
      if(Nyx::theGhostFDMPC())
      	Nyx::theGhostFDMPC()->DepositFDMParticles(fdmreal,fdmimag,lev,a_new);
      if(Nyx::theVirtFDMPC())
      	Nyx::theVirtFDMPC()->DepositFDMParticles(fdmreal,fdmimag,lev,a_new);
      if(Nyx::theFDMwkbPC())
      	Nyx::theFDMwkbPC()->DepositFDMParticles(fdmreal,fdmimag,lev,a_new,Nyx::theta_fdm,hbaroverm);
      if(Nyx::theGhostFDMwkbPC())
      	Nyx::theGhostFDMwkbPC()->DepositFDMParticles(fdmreal,fdmimag,lev,a_new,Nyx::theta_fdm,hbaroverm);
      if(Nyx::theVirtFDMwkbPC())
      	Nyx::theVirtFDMwkbPC()->DepositFDMParticles(fdmreal,fdmimag,lev,a_new,Nyx::theta_fdm,hbaroverm);
      if(Nyx::theFDMphasePC())
      	Nyx::theFDMphasePC()->DepositFDMParticles(fdmreal,fdmimag,lev,a_new,Nyx::theta_fdm,hbaroverm);
      if(Nyx::theGhostFDMphasePC())
      	Nyx::theGhostFDMphasePC()->DepositFDMParticles(fdmreal,fdmimag,lev,a_new,Nyx::theta_fdm,hbaroverm);
      if(Nyx::theVirtFDMphasePC())
      	Nyx::theVirtFDMphasePC()->DepositFDMParticles(fdmreal,fdmimag,lev,a_new,Nyx::theta_fdm,hbaroverm);

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
          BL_FORT_PROC_CALL(FORT_FDM_FIELDS, fort_fdm_fields)
            (BL_TO_FORTRAN(Ax_new[fpi]));
          if (Ax_new[fpi].contains_nan())
	    amrex::Abort("Nans in state just after FDM density update");
        }
      } else if (levelmethod[lev]==CWlevel){

      amrex::Print() << "levelSteps " << lev << " "<< parent->levelSteps(lev) << '\n';

      //Define neccessary number of ghost cells                                                                                                                                                                 
      //int ng = parent->nCycle(lev)+2.0*ceil(Nyx::sigma_fdm*Nyx::theta_fdm/get_level(lev).Geom().CellSize()[0]);
      int ng = parent->nCycle(lev) + 2.0*ceil(Nyx::theta_fdm);

      //Initialize MultiFabs                                                                                                                                                                                    
      MultiFab& Ax_new = get_level(lev).get_new_data(Axion_Type);
      Ax_new.setVal(0.);
      MultiFab fdmreal(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
      fdmreal.setVal(0.);
      MultiFab fdmimag(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
      fdmimag.setVal(0.);
      MultiFab fdmdens(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng); // multifab for N-body density
      fdmdens.setVal(0.); 

      //Deposit CWA                                                                                                                                         
      if(Nyx::theFDMphasePC())
      	Nyx::theFDMphasePC()->DepositFDMParticlesCWA(fdmreal,fdmimag,lev,a_new,Nyx::theta_fdm,hbaroverm);
      if(Nyx::theGhostFDMphasePC())
      	Nyx::theGhostFDMphasePC()->DepositFDMParticlesCWA(fdmreal,fdmimag,lev,a_new,Nyx::theta_fdm,hbaroverm);
      if(Nyx::theVirtFDMphasePC())
      	Nyx::theVirtFDMphasePC()->DepositFDMParticlesCWA(fdmreal,fdmimag,lev,a_new,Nyx::theta_fdm,hbaroverm);

      // Assign N-body density and add it to Ax_new

      // if(Nyx::theFDMphasePC()){
      // 	Nyx::theFDMphasePC()->AssignDensitySingleLevel(fdmdens,lev); 
      // 	MultiFab::Add(Ax_new,fdmdens, 0, Nyx::AxDens, 1, 0);
      // }
      // if(Nyx::theGhostFDMphasePC()){
      // 	Nyx::theGhostFDMphasePC()->AssignDensitySingleLevel(fdmdens,lev);
      // 	MultiFab::Add(Ax_new,fdmdens, 0, Nyx::AxDens, 1, 0);
      // }
      // if(Nyx::theVirtFDMphasePC()){
      // 	Nyx::theVirtFDMphasePC()->AssignDensitySingleLevel(fdmdens,lev);
      // 	MultiFab::Add(Ax_new,fdmdens, 0, Nyx::AxDens, 1, 0);
      // }

      // Ax_new.FillBoundary(Nyx::AxDens, 1, parent->Geom(lev).periodicity());


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
          BL_FORT_PROC_CALL(FORT_FDM_FIELDS, fort_fdm_fields) // fort_fdm_fields2 if N-body is used 
            (BL_TO_FORTRAN(Ax_new[fpi]));
          if (Ax_new[fpi].contains_nan())
	    amrex::Abort("Nans in state just after FDM density update");
        }

      } else
	continue;
    }
#endif

    // Redistribution happens in post_timestep
    return dt;
}
#endif
