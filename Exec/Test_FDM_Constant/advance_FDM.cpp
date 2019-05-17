#ifdef FDM

#include "Nyx.H"
#include "Nyx_F.H"
#include <AMReX_Particles_F.H>
#include <AMReX_ParallelDescriptor.H>
using namespace amrex;


#ifdef GRAVITY
#	include "Gravity.H"
#	include <Gravity_F.H>
#endif

using std::string;

amrex::Real
Nyx::advance_FDM (amrex::Real time,
                              amrex::Real dt,
                              int  iteration,
                              int  ncycle)
{
    bool show_timings=false;
    //I put this here; its value is set differently in Nyx_advance (advance_hydro_plus_particles), but in the old version of the code, it was set to 3.
    int grav_n_grow = 3;
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
        //setup_virtual_particles();
        //
        // Setup ghost particles for use in finer levels. Note that Ghost particles
        // that will be used by this level have already been created, the
        // particles being set here are only used by finer levels.
        //
       // for(int lev = level; lev <= finest_level_to_advance && lev < finest_level; lev++)
       // {
       //    get_level(lev).setup_ghost_particles();
       // }
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

            MultiFab& S_old = get_level(lev).get_old_data(State_Type);
            MultiFab& S_new = get_level(lev).get_new_data(State_Type);
            MultiFab::Copy(S_new, S_old, 0, 0, S_old.nComp(), 0);

            MultiFab& D_old = get_level(lev).get_old_data(DiagEOS_Type);
            MultiFab& D_new = get_level(lev).get_new_data(DiagEOS_Type);
            MultiFab::Copy(D_new, D_old, 0, 0, D_old.nComp(), 0);
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
#ifdef CGRAV
        MultiFab& phi_old = get_level(level).get_old_data(PhiGrav_Type);
        MultiFab& phi_new = get_level(level).get_new_data(PhiGrav_Type);

        prescribe_grav_potential(phi_old, Geom() , level, finest_level);

        MultiFab::Copy(phi_new, phi_old, 0, 0, phi_old.nComp(), 0);


#else
        gravity->multilevel_solve_for_old_phi(level, finest_level,
                                              use_previous_phi_as_guess);
#endif
    }
//    //
//    // Advance Particles
//    //
//    if (Nyx::theActiveParticles().size() > 0)
//    {
//        // Advance the particle velocities to the half-time and the positions to the new time
//        // We use the cell-centered gravity to correctly interpolate onto particle locations
//        if (particle_move_type == "Gravitational")
//        {
//            const amrex::Real a_half = 0.5 * (a_old + a_new);
//
//            if (particle_verbose && ParallelDescriptor::IOProcessor())
//                std::cout << "moveKickDrift ... updating particle positions and velocity\n";
//
//            for (int lev = level; lev <= finest_level_to_advance; lev++)
//            {
//                // We need grav_n_grow grow cells to track boundary particles
//                const BoxArray& ba = get_level(lev).get_new_data(State_Type).boxArray();
//                MultiFab grav_vec_old(ba, BL_SPACEDIM, grav_n_grow);
//                get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old, time);
//
//                for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
//                    Nyx::theActiveParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half);
//
//                // Only need the coarsest virtual particles here.
//                if (lev == level && level < finest_level)
//                    for (int i = 0; i < Nyx::theVirtualParticles().size(); i++)
//                        Nyx::theVirtualParticles()[i]->moveKickDrift(grav_vec_old, level, dt, a_old, a_half);
//
//                // Miiiight need all Ghosts
//                for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
//                    Nyx::theGhostParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_new, a_half);
//            }
//        }
//    }
//
#endif
//TODO_JENS: here is the spot to put the hooks for calling the FFT solver.
    //Advance Axions

    for (int lev = level; lev <= finest_level_to_advance; lev++)
        if(lev==0)
            //here is the hook:
            get_level(lev).advance_FDM_FFT(time, dt, a_old, a_new);
        else
            get_level(lev).advance_FDM_FD(time, dt, a_old, a_new);

    // Always average down from finer to coarser.
    for (int lev = finest_level_to_advance-1; lev >= level; lev--)
        get_level(lev).average_down();

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
//TODO: same here
//
//    // Must average down again after doing the gravity correction;
//    //      always average down from finer to coarser.
//    //    for (int lev = finest_level_to_advance-1; lev >= level; lev--)
//    //        get_level(lev).average_down();
//
//    if (Nyx::theActiveParticles().size() > 0)
//    {
//        // Advance the particle velocities by dt/2 to the new time. We use the
//        // cell-centered gravity to correctly interpolate onto particle
//        // locations.
//        if (particle_move_type == "Gravitational")
//        {
//            const amrex::Real a_half = 0.5 * (a_old + a_new);
//
//            if (particle_verbose && ParallelDescriptor::IOProcessor())
//                std::cout << "moveKick ... updating velocity only\n";
//
//            for (int lev = level; lev <= finest_level_to_advance; lev++)
//            {
//                const BoxArray& ba = get_level(lev).get_new_data(State_Type).boxArray();
//                MultiFab grav_vec_new(ba, BL_SPACEDIM, grav_n_grow);
//                get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new, cur_time);
//
//                for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
//                    Nyx::theActiveParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);
//
//                // Virtual particles will be recreated, so we need not kick them.
//
//                // Ghost particles need to be kicked except during the final iteration.
//                if (iteration != ncycle)
//                    for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
//                        Nyx::theGhostParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);
//            }
//        }
//    }

    // //we need to get new grids, since it implicitely fills the new part of the
    // //Gravity_Type, which we need for the interpolated values of G at higher
    // //levels when filling ghosts from lower levels.
    // const auto& dm = get_level(level).get_new_data(State_Type).DistributionMap();
    // MultiFab grav_vec_new(grids, dm, BL_SPACEDIM, grav_n_grow);
    // get_level(level).gravity->get_new_grav_vector(level, grav_vec_new, cur_time);

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Time  after moveKick " << end << '\n';
    }
#endif

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
