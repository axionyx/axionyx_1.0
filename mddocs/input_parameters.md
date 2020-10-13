## Input Parameters
This section briefly describes the parameter available in the inputs file.

`nyx.v`- _(optional)_ this is an example of a docstring for an input parameter. This line lives in Source/Nyx.cpp, line 258  
`stopwatch.level`- _(optional)_ level of verbosity for Stopwatch.  

####Resolution parameters 

`amr.n_cell           = 64 64 64`  - the number of grid points along one dimension, for a cubic 3D grid this is repeated 3-times. 
`amr.max_grid_size    = 16` - the maximum size of a subgrid. 
`amr.blocking_factor  = 16` - the minumum size of a subgrid when using refinement. Without refinement, this can be set to the same value as `amr.max_grid_size`.  

There are some restrictions to what these parameters are allowed to be.

Rules: 
1) number of MPI ranks (m) has to be 2^(3*i), where i=0,1,2,3,... so m = 1,8,64,512,...
2) amr.max_grid_size should always be set to (N^3 /m)^1/3, where N = amr.n_cell

Explanation: 
The 3-dimensional grid with the size N^3 is distributed to m MPI processors. But pencil decomposition that is used to parallelise the Fourier transform also requires that the grid be split into powers of 2 along one dimension. So in 3 dimensions, m has to be of the form 2^(3*i).

It follows then that `amr.max_grid_size` or the size of the subgrid along one dimension is (total volume / number of MPI ranks)^1/3. 
 
So in the above example with 

`amr.n_cell = 64 64 64` 

one woudl use 
`amr.max_grid_size = 64` for 1 MPI rank
`amr.max_grid_size = 32` for  8 ranks
`amr.max_grid_size = 16` for 64 ranks and so on. 


### Undocumented Input Parameters

`nyx.init_shrink`-   
`nyx.cfl`-   
`nyx.change_max`-   
`nyx.fixed_dt`-   
`nyx.initial_dt`-   
`nyx.sum_interval`-   
`nyx.do_reflux`-   
`nyx.dt_cutoff`-   
`nyx.vonNeumann_dt`-   
`nyx.m_tt`-   
`nyx.theta_fdm`-   
`nyx.sigma_fdm`-   
`nyx.alpha_fdm`-   
`nyx.wkb_approx`-   
`nyx.order`-   
`nyx.beam_cfl`-   
`particles.part_size`-   
`amr.n_cell`-   
`nyx.initial_z`-   
`nyx.final_a`-   
`nyx.final_z`-   
`nyx.relative_max_change_a`-   
`nyx.absolute_max_change_a`-   
`nyx.dt_binpow`-   
`nyx.plot_rank`-   
`nyx.plot_X`-   
`nyx.lo_bc`-   
`nyx.hi_bc`-   
`max_step`-   
`strt_time`-   
`stop_time`-   
`how`-   
`nyx.do_dm_particles`-   
`nyx.particle_init_type`-   
`nyx.particle_move_type`-   
`nyx.init_with_sph_particles`-   
`nyx.particle_initrandom_serialize`-   
`nyx.particle_initrandom_count`-   
`nyx.particle_initrandom_count_per_box`-   
`nyx.particle_initrandom_mass`-   
`nyx.particle_initrandom_iseed`-   
`nyx.particle_skip_factor`-   
`nyx.ascii_particle_file`-   
`nyx.sph_particle_file`-   
`nyx.binary_particle_file`-   
`nyx.agn_particle_file`-   
`nyx.neutrino_particle_file`-   
`nyx.num_particle_fdm`-   
`nyx.num_particle_dm`-   
`nyx.write_particle_density_at_init`-   
`nyx.write_coarsened_particles`-   
`particles.v`-   
`particles.replicate`-   
`particles.cfl`-   
`particles.neutrino_cfl`-   
`particles.dm_particle_output_file`-   
`particles.agn_particle_output_file`-   
`particles.fdm_particle_output_file`-   
`forcing.v`-   
`forcing.seed`-   
`forcing.profile`-   
`forcing.soln_weight`-   
`forcing.alpha`-   
`forcing.band_width`-   
`forcing.intgr_vel`-   
`forcing.auto_corrl`-   
`gravity.gravity_type`-   
`gravity.v`-   
`gravity.no_sync`-   
`gravity.no_composite`-   
`gravity.dirichlet_bcs`-   
`gravity.mlmg_max_fmg_iter`-   
`gravity.mlmg_agglomeration`-   
`gravity.mlmg_consolidation`-   
`gravity.ml_tol`-   
`gravity.sl_tol`-   
`gravity.delta_tol`-   
`cosmo.initDirName`-   
`cosmo.nuInitDirName`-   
`nyx.n_particles`-   
`nyx.do_santa_barbara`-   
`nyx.init_sb_vels`-   
`nyx.do_readin_ics`-   
`nyx.readin_ics_fname`-   
`nyx.eos_nr_eps`-   
`nyx.vode_rtol`-   
`nyx.vode_atol_scaled`-   
