#include <iomanip>
#include <Nyx.H>

#ifdef GRAVITY
#include <Gravity.H>
#include <Gravity_F.H> //needed for get_grav_constant, but there might be a better way
#endif

#include <Nyx_F.H>
#include <AMReX_Particles_F.H>

#include "fdm_F.H"

using namespace amrex;

namespace
{
    bool virtual_particles_set = false;
    
    std::string ascii_particle_file;
    std::string binary_particle_file;
    std::string    sph_particle_file;

#ifdef AGN
    std::string agn_particle_file;
#endif

#ifdef NEUTRINO_PARTICLES
    std::string neutrino_particle_file;
#endif

#ifdef FDM
    std::string fdm_particle_file;
#endif

  // const std::string chk_particle_file("DM");
    const std::string dm_chk_particle_file("DM");
    const std::string agn_chk_particle_file("AGN");
#ifdef FDM
    const std::string fdm_chk_particle_file("FDM");
#endif
    //
    // We want to call this routine when on exit to clean up particles.
    //

    //
    // Array of containers for all active particles
    //
    Vector<NyxParticleContainerBase*> ActiveParticles;
    //
    // Array of containers for all virtual particles
    //
    Vector<NyxParticleContainerBase*> VirtualParticles;
    //
    // Array of containers for all ghost particles
    //
    Vector<NyxParticleContainerBase*> GhostParticles;

    //
    // Containers for the real "active" Particles
    //
    DarkMatterParticleContainer* DMPC = 0;
    StellarParticleContainer*     SPC = 0;
#ifdef AGN
    AGNParticleContainer*         APC = 0;
#endif
#ifdef NEUTRINO_PARTICLES
    NeutrinoParticleContainer*    NPC = 0;
#endif
#ifdef FDM
    FDMwkbParticleContainer* FDMwkbPC = 0;
    FDMParticleContainer*       FDMPC = 0;
#endif
    //
    // This is only used as a temporary container for
    //  reading in SPH particles and using them to
    //  initialize the density and velocity field on the
    //  grids.
    //
    DarkMatterParticleContainer*        SPHPC = 0;
    //
    // Container for temporary, virtual Particles
    //
    DarkMatterParticleContainer* VirtPC  = 0;
    StellarParticleContainer*    VirtSPC = 0;
#ifdef AGN
    AGNParticleContainer*        VirtAPC = 0;
#endif
#ifdef NEUTRINO_PARTICLES
    NeutrinoParticleContainer*   VirtNPC = 0;
#endif
#ifdef FDM
    FDMwkbParticleContainer* VirtFDMwkbPC = 0;
    FDMParticleContainer*       VirtFDMPC = 0;
#endif
    //
    // Container for temporary, ghost Particles
    //
    DarkMatterParticleContainer* GhostPC  = 0;
    StellarParticleContainer*    GhostSPC = 0;
#ifdef AGN
    AGNParticleContainer*        GhostAPC = 0;
#endif
#ifdef NEUTRINO_PARTICLES
    NeutrinoParticleContainer*   GhostNPC = 0;
#endif
#ifdef FDM
    FDMwkbParticleContainer* GhostFDMwkbPC = 0;
    FDMParticleContainer*       GhostFDMPC = 0;
#endif

    void RemoveParticlesOnExit ()
    {
        for (int i = 0; i < ActiveParticles.size(); i++)
        {
            delete ActiveParticles[i];
            ActiveParticles[i] = 0;
        }
        for (int i = 0; i < GhostParticles.size(); i++)
        {
            delete GhostParticles[i];
            GhostParticles[i] = 0;
        }
        for (int i = 0; i < VirtualParticles.size(); i++)
        {
            delete VirtualParticles[i];
            VirtualParticles[i] = 0;
        }
#ifdef FDM
	if(FDMwkbPC)
	  delete FDMwkbPC;
	if(GhostFDMwkbPC)
	  delete GhostFDMwkbPC;
	if(VirtFDMwkbPC)
	  delete VirtFDMwkbPC;
	if(FDMPC)
	  delete FDMPC;
	if(GhostFDMPC)
	  delete GhostFDMPC;
	if(VirtFDMPC)
	  delete VirtFDMPC;
#endif
    }
}

bool Nyx::do_dm_particles = false;
int Nyx::num_particle_ghosts = 1;
int Nyx::particle_skip_factor = 1;

std::string Nyx::particle_init_type = "";
std::string Nyx::particle_move_type = "";

// Allows us to output particles in the plotfile
//   in either single (IEEE32) or double (NATIVE) precision.  
// Particles are always written in double precision
//   in the checkpoint files.

bool Nyx::particle_initrandom_serialize = false;
Real Nyx::particle_initrandom_mass;
long Nyx::particle_initrandom_count;
long Nyx::particle_initrandom_count_per_box;
int  Nyx::particle_initrandom_iseed;

int Nyx::particle_verbose               = 1;
int Nyx::write_particle_density_at_init = 0;
int Nyx::write_coarsened_particles      = 0;
Real Nyx::particle_cfl = 0.5;
#ifdef NEUTRINO_PARTICLES
Real Nyx::neutrino_cfl = 0.5;
#endif
#ifdef FDM
long Nyx::num_particle_fdm;
long Nyx::num_particle_dm;
#endif

IntVect Nyx::Nrep;

Vector<NyxParticleContainerBase*>&
Nyx::theActiveParticles ()
{
    return ActiveParticles;
}

Vector<NyxParticleContainerBase*>&
Nyx::theGhostParticles ()
{
    return GhostParticles;
}

Vector<NyxParticleContainerBase*>&
Nyx::theVirtualParticles ()
{
    return VirtualParticles;
}

DarkMatterParticleContainer*
Nyx::theDMPC ()
{
    return DMPC;
}

DarkMatterParticleContainer*
Nyx::theVirtPC ()
{
    return VirtPC;
}
DarkMatterParticleContainer*
Nyx::theGhostPC ()
{
    return GhostPC;
}

StellarParticleContainer* 
Nyx::theSPC ()
{
      return SPC;
}
StellarParticleContainer* 
Nyx::theVirtSPC ()
{
      return VirtSPC;
}
StellarParticleContainer* 
Nyx::theGhostSPC ()
{
      return GhostSPC;
}

#ifdef AGN
AGNParticleContainer* 
Nyx::theAPC ()
{
      return APC;
}
AGNParticleContainer* 
Nyx::theVirtAPC ()
{
      return VirtAPC;
}
AGNParticleContainer* 
Nyx::theGhostAPC ()
{
      return GhostAPC;
}
#endif

#ifdef NEUTRINO_PARTICLES
NeutrinoParticleContainer* 
Nyx::theNPC ()
{
      return NPC;
}
NeutrinoParticleContainer* 
Nyx::theVirtNPC ()
{
      return VirtNPC;
}
NeutrinoParticleContainer* 
Nyx::theGhostNPC ()
{
      return GhostNPC;
}
#endif

#ifdef FDM
FDMParticleContainer*
Nyx::theFDMPC ()
{
  return FDMPC;
}
FDMParticleContainer*
Nyx::theVirtFDMPC ()
{
  return VirtFDMPC;
}
FDMParticleContainer*
Nyx::theGhostFDMPC ()
{
  return GhostFDMPC;
}
FDMwkbParticleContainer*
Nyx::theFDMwkbPC ()
{
  return FDMwkbPC;
}
FDMwkbParticleContainer*
Nyx::theVirtFDMwkbPC ()
{
  return VirtFDMwkbPC;
}
FDMwkbParticleContainer*
Nyx::theGhostFDMwkbPC ()
{
  return GhostFDMwkbPC;
}
#endif

void
Nyx::read_particle_params ()
{

    ParmParse pp("nyx");
    pp.query("do_dm_particles", do_dm_particles);
#ifdef FDM
    if(partlevel)
      //If Gaussian Beam Method is chosen for FDM evolution, we need dm_particles for the construction of the gravitational potential.
      do_dm_particles = 1;
#endif
#ifdef AGN
    pp.get("particle_init_type", particle_init_type);
    pp.get("particle_move_type", particle_move_type);
#else
    if (do_dm_particles)
    {
        pp.get("particle_init_type", particle_init_type);
        pp.get("particle_move_type", particle_move_type);
        pp.query("init_with_sph_particles", init_with_sph_particles);
    }
#endif

#ifdef GRAVITY
    if (!do_grav && particle_move_type == "Gravitational")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR:: doesnt make sense to have do_grav=false but move_type = Gravitational" << std::endl;
        amrex::Error();
    }
#endif

    pp.query("particle_initrandom_serialize", particle_initrandom_serialize);
    pp.query("particle_initrandom_count", particle_initrandom_count);
    pp.query("particle_initrandom_count_per_box", particle_initrandom_count_per_box);
    pp.query("particle_initrandom_mass", particle_initrandom_mass);
    pp.query("particle_initrandom_iseed", particle_initrandom_iseed);
    pp.query("particle_skip_factor", particle_skip_factor);
    pp.query("ascii_particle_file", ascii_particle_file);

    // Input error check
    if (do_dm_particles && !ascii_particle_file.empty() && particle_init_type != "AsciiFile")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not AsciiFile but you specified ascii_particle_file" << std::endl;;
        amrex::Error();
    }

    pp.query("sph_particle_file", sph_particle_file);

    // Input error check
    if (init_with_sph_particles != 1 && !sph_particle_file.empty())
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::init_with_sph_particles is not 1 but you specified sph_particle_file" << std::endl;;
        amrex::Error();
    }

    // Input error check
    if (init_with_sph_particles == 1 && sph_particle_file.empty())
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::init_with_sph_particles is 1 but you did not specify sph_particle_file" << std::endl;;
        amrex::Error();
    }

    pp.query("binary_particle_file", binary_particle_file);

    // Input error check
    if (!binary_particle_file.empty() && (particle_init_type != "BinaryFile" &&
                                          particle_init_type != "BinaryMetaFile" && 
					  particle_init_type != "BinaryMortonFile"))
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not BinaryFile, BinaryMetaFile, or BinaryMortonFile but you specified binary_particle_file" << std::endl;
        amrex::Error();
    }

#ifdef AGN
    pp.query("agn_particle_file", agn_particle_file);
    if (!agn_particle_file.empty() && particle_init_type != "AsciiFile")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not AsciiFile but you specified agn_particle_file" << std::endl;;
        amrex::Error();
    }
#endif

#ifdef NEUTRINO_PARTICLES
    pp.query("neutrino_particle_file", neutrino_particle_file);
    if (!neutrino_particle_file.empty() && particle_init_type != "AsciiFile")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not AsciiFile but you specified neutrino_particle_file" << std::endl;;
        amrex::Error();
    }
#endif

#ifdef FDM
    if(partlevel){
      pp.get("num_particle_fdm", num_particle_fdm);
      pp.get("num_particle_dm", num_particle_dm);
    }
#endif

    pp.query("write_particle_density_at_init", write_particle_density_at_init);
    pp.query("write_coarsened_particles", write_coarsened_particles);
    //
    // Control the verbosity of the Particle class
    //
    ParmParse ppp("particles");
    ppp.query("v", particle_verbose);

    for (int i = 0; i < BL_SPACEDIM; i++) Nrep[i] = 1; // Initialize to one (no replication)
    ppp.query("replicate",Nrep);
    //
    // Set the cfl for particle motion (fraction of cell that a particle can
    // move in a timestep).
    //
    ppp.query("cfl", particle_cfl);
#ifdef NEUTRINO_PARTICLES
    ppp.query("neutrino_cfl", neutrino_cfl);
#endif
}

void
Nyx::init_particles ()
{
    BL_PROFILE("Nyx::init_particles()");

    if (level > 0)
        return;

#ifdef FDM
    const Real a = get_comoving_a(state[State_Type].curTime());
#endif
    //
    // Need to initialize particles before defining gravity.
    //
    if (do_dm_particles)
    {
        BL_ASSERT (DMPC == 0);
        DMPC = new DarkMatterParticleContainer(parent);
        ActiveParticles.push_back(DMPC); 

        if (init_with_sph_particles == 1)
           SPHPC = new DarkMatterParticleContainer(parent);

	if (parent->subCycle())
	{
	    VirtPC = new DarkMatterParticleContainer(parent);
            VirtualParticles.push_back(VirtPC); 

	    GhostPC = new DarkMatterParticleContainer(parent);
            GhostParticles.push_back(GhostPC); }

        //
        // Make sure to call RemoveParticlesOnExit() on exit.
        //
        amrex::ExecOnFinalize(RemoveParticlesOnExit);
        //
        // 2 gives more stuff than 1.
        //
        DMPC->SetVerbose(particle_verbose);

        DarkMatterParticleContainer::ParticleInitData pdata = {particle_initrandom_mass};

        if (particle_init_type == "Random")
        {
            if (particle_initrandom_count <= 0)
            {
                amrex::Abort("Nyx::init_particles(): particle_initrandom_count must be > 0");
            }
            if (particle_initrandom_iseed <= 0)
            {
                amrex::Abort("Nyx::init_particles(): particle_initrandom_iseed must be > 0");
            }

            if (verbose)
            {
                amrex::Print() << "\nInitializing DM with cloud of "
                               << particle_initrandom_count
                               << " random particles with initial seed: "
                               << particle_initrandom_iseed << "\n\n";
            }

            DMPC->InitRandom(particle_initrandom_count,
                             particle_initrandom_iseed, pdata,
                             particle_initrandom_serialize);
        }
        else if (particle_init_type == "RandomPerBox")
        {
            if (particle_initrandom_count_per_box <= 0)
            {
                amrex::Abort("Nyx::init_particles(): particle_initrandom_count_per_box must be > 0");
            }
            if (particle_initrandom_iseed <= 0)
            {
                amrex::Abort("Nyx::init_particles(): particle_initrandom_iseed must be > 0");
            }

            if (verbose)
                amrex::Print() << "\nInitializing DM with of " << particle_initrandom_count_per_box
                               << " random particles per box with initial seed: "
                               << particle_initrandom_iseed << "\n\n";

            DMPC->InitRandomPerBox(particle_initrandom_count_per_box,
                                   particle_initrandom_iseed, pdata);

        }
        else if (particle_init_type == "RandomPerCell")
        {
            if (verbose)
                amrex::Print() << "\nInitializing DM with 1 random particle per cell " << "\n";

            int n_per_cell = 1;
            DMPC->InitNRandomPerCell(n_per_cell, pdata);
        }
        else if (particle_init_type == "AsciiFile")
        {
            if (verbose)
            {
                amrex::Print() << "\nInitializing DM particles from \""
                               << ascii_particle_file << "\" ...\n\n";
                if (init_with_sph_particles == 1)
                    amrex::Print() << "\nInitializing SPH particles from ascii \""
                                   << sph_particle_file << "\" ...\n\n";
            }
            //
            // The second argument is how many Reals we read into `m_data[]`
            // after reading in `m_pos[]`. Here we're reading in the particle
            // mass and velocity.
            //
            DMPC->InitFromAsciiFile(ascii_particle_file, BL_SPACEDIM + 1, &Nrep);
            if (init_with_sph_particles == 1)
                SPHPC->InitFromAsciiFile(ascii_particle_file, BL_SPACEDIM + 1, &Nrep);
        }
        else if (particle_init_type == "BinaryFile")
        {
            if (verbose)
            {
                amrex::Print() << "\nInitializing DM particles from \""
                               << binary_particle_file << "\" ...\n\n";
                if (init_with_sph_particles == 1)
                    amrex::Print() << "\nInitializing SPH particles from binary \""
                                   << sph_particle_file << "\" ...\n\n";
            }
            //
            // The second argument is how many Reals we read into `m_data[]`
            // after reading in `m_pos[]`. Here we're reading in the particle
            // mass and velocity.
            //
            DMPC->InitFromBinaryFile(binary_particle_file, BL_SPACEDIM + 1);
            if (init_with_sph_particles == 1)
                SPHPC->InitFromBinaryFile(ascii_particle_file, BL_SPACEDIM + 1);
        }
        else if (particle_init_type == "BinaryMetaFile")
        {
            if (verbose)
            {
                amrex::Print() << "\nInitializing DM particles from meta file\""
                               << binary_particle_file << "\" ...\n\n";
                if (init_with_sph_particles == 1)
                    amrex::Print() << "\nInitializing SPH particles from meta file\""
                                   << sph_particle_file << "\" ...\n\n";
            }
            //
            // The second argument is how many Reals we read into `m_data[]`
            // after reading in `m_pos[]` in each of the binary particle files.
            // Here we're reading in the particle mass and velocity.
            //
            DMPC->InitFromBinaryMetaFile(binary_particle_file, BL_SPACEDIM + 1);
            if (init_with_sph_particles == 1)
                SPHPC->InitFromBinaryMetaFile(sph_particle_file, BL_SPACEDIM + 1);
        }
        else if (particle_init_type == "BinaryMortonFile")
        {
	  if (verbose)
            {
	      amrex::Print() << "\nInitializing DM particles from morton-ordered binary file\""
			     << binary_particle_file << "\" ...\n\n";
	      if (init_with_sph_particles == 1)
		amrex::Error("Morton-ordered input is not supported for sph particles.");
            }
            //
            // The second argument is how many Reals we read into `m_data[]`
            // after reading in `m_pos[]` in each of the binary particle files.
            // Here we're reading in the particle mass and velocity.
            //
	  DMPC->InitFromBinaryMortonFile(binary_particle_file,
					 BL_SPACEDIM + 1,
					 particle_skip_factor);
        }
#ifdef FDM
	else if (particle_init_type == "GaussianBeams" && partlevel)
	  {
	    if (verbose)
	      {
		amrex::Print() << "\nInitializing DM particles for FDM gravitational potential...\n\n";
		if (init_with_sph_particles == 1)
		  amrex::Error("FDM computations are not supported for sph particles.");
	      }
	    if(num_particle_dm > 0)
	      DMPC->InitGaussianBeams(num_particle_dm, level, parent->initialBaLevels()+1, meandens, alpha_fdm, a);
	    else
	      amrex::Error("\nNeed num_particle_dm > 0 for DM InitGaussianBeams!\n\n");
	  }
#endif
        else
        {
            amrex::Error("not a valid input for nyx.particle_init_type");
        }

        if (write_particle_density_at_init == 1)
        {
            MultiFab particle_mf(grids,dmap,1,1);
            DMPC->AssignDensitySingleLevel(particle_mf,0,1,0);

            writeMultiFabAsPlotFile("ParticleDensity", particle_mf, "density");
            exit(0);
        }

        if (write_coarsened_particles)
        {
            DMPC->WriteCoarsenedAsciiFile("coarse_particle_file.ascii");
            exit(0);
        }
    }
#ifdef AGN
    {
        // Note that we don't initialize any actual AGN particles here, we just create the container.
        BL_ASSERT (APC == 0);
        APC = new AGNParticleContainer(parent, num_particle_ghosts);
        ActiveParticles.push_back(APC); 

	if (parent->subCycle())
	{
	    VirtAPC = new AGNParticleContainer(parent, num_particle_ghosts);
            VirtualParticles.push_back(VirtAPC); 

	    GhostAPC = new AGNParticleContainer(parent, num_particle_ghosts);
            GhostParticles.push_back(GhostAPC); 
	}
        //
        // 2 gives more stuff than 1.
        //
        APC->SetVerbose(particle_verbose);
    }
#endif
#ifdef NEUTRINO_PARTICLES
    {
        BL_ASSERT (NPC == 0);
        NPC = new NeutrinoParticleContainer(parent);
        ActiveParticles.push_back(NPC); 

        // Set the relativistic flag to 1 so that the density will be gamma-weighted
        //     when returned by AssignDensity calls.
        NPC->SetRelativistic(1);

        // We must set the value for csquared which is used in computing gamma = 1 / sqrt(1-vsq/csq)
        // Obviously this value is just a place-holder for now.
        NPC->SetCSquared(1.);

	if (parent->subCycle())
	{
	    VirtNPC = new NeutrinoParticleContainer(parent);
            VirtualParticles.push_back(VirtNPC); 

            // Set the relativistic flag to 1 so that the density will be gamma-weighted
            //     when returned by AssignDensity calls.
            VirtNPC->SetRelativistic(1);

	    GhostNPC = new NeutrinoParticleContainer(parent);
            GhostParticles.push_back(GhostNPC); 

            // Set the relativistic flag to 1 so that the density will be gamma-weighted
            //     when returned by AssignDensity calls.
            GhostNPC->SetRelativistic(1);
	}
        //
        // Make sure to call RemoveParticlesOnExit() on exit.
        //   (if do_dm_particles then we have already called ExecOnFinalize)
        //
        if (!do_dm_particles)
            amrex::ExecOnFinalize(RemoveParticlesOnExit);
        //
        // 2 gives more stuff than 1.
        //
        NPC->SetVerbose(particle_verbose);
        if (particle_init_type == "AsciiFile")
        {
            if (verbose)
                amrex::Print() << "\nInitializing Neutrino particles from \""
                               << neutrino_particle_file << "\" ...\n\n";
            //
            // The second argument is how many Reals we read into `m_data[]`
            // after reading in `m_pos[]`. Here we're reading in the particle
            // mass, velocity and angles.
            //
            NPC->InitFromAsciiFile(neutrino_particle_file, 2*BL_SPACEDIM + 1, &Nrep);
        }
        else if (particle_init_type == "BinaryFile")
        {
            if (verbose)
            {
                amrex::Print() << "\nInitializing Neutrino particles from \""
                               << neutrino_particle_file << "\" ...\n\n";
            }
            //
            // The second argument is how many Reals we read into `m_data[]`
            // after reading in `m_pos[]`. Here we're reading in the particle
            // mass and velocity.
            //
            NPC->InitFromBinaryFile(neutrino_particle_file, BL_SPACEDIM + 1);
        }

        else
        {
            amrex::Error("for right now we only init Neutrino particles with ascii or binary");
        }
    }
#endif
#ifdef FDM
    if(partlevel){
    // if(false){
    if(wkb_approx)
    {
      BL_ASSERT (FDMwkbPC == 0);
      FDMwkbPC = new FDMwkbParticleContainer(parent);

      if (parent->subCycle())
        {
      	  VirtFDMwkbPC = new FDMwkbParticleContainer(parent);

      	  GhostFDMwkbPC = new FDMwkbParticleContainer(parent);
        }

      //                                                                                                                                                                                                         
      // Make sure to call RemoveParticlesOnExit() on exit.                                                                                                                                                      
      //   (if do_dm_particles then we have already called ExecOnFinalize)                                                                                                                                       
      //                                                                                                                                                                                                         
      if (!do_dm_particles)
	amrex::ExecOnFinalize(RemoveParticlesOnExit);
      //                                                                                                                                                                                                         
      // 2 gives more stuff than 1.                                                                                                                                                                              
      //                                                                                                                                                                                                         
      FDMwkbPC->SetVerbose(particle_verbose);

      if (particle_init_type == "AsciiFile")
        {
	  if (verbose && ParallelDescriptor::IOProcessor())
            {
	      amrex::Print() << "\nInitializing FDM particles from \""
			     << ascii_particle_file << "\" ...\n\n";
            }
	  //                                                                                                                                                                                                     
	  // The second argument is how many Reals we read into `m_data[]`                                                                                                                                       
	  // after reading in `m_pos[]`. Here we're reading in the particle                                                                                                                                      
	  // mass and velocity. This might have to be extended to all Gauss 
	  // Beam arguments.                                                                                                                                                               
	  FDMwkbPC->InitFromAsciiFile(ascii_particle_file, BL_SPACEDIM + 1, &Nrep);
        }
      else if (particle_init_type == "BinaryFile")
        {
	  if (verbose)
            {
	      amrex::Print() << "\nInitializing FDM particles from \""
			     << binary_particle_file << "\" ...\n\n";
            }
	  //                                                                                                                                                                                                     
	  // The second argument is how many Reals we read into `m_data[]`                                                                                                                                       
	  // after reading in `m_pos[]`. Here we're reading in the particle                                                                                                                                      
	  // mass and velocity. This might have to be extended to all Gauss                                                                                                                          
          // Beam arguments.                                                                                                                                                                           
	  FDMwkbPC->InitFromBinaryFile(binary_particle_file, BL_SPACEDIM + 1);
        }
      else if (particle_init_type == "BinaryMetaFile")
        {
	  if (verbose)
            {
	      amrex::Print() << "\nInitializing FDM particles from meta file\""
			     << binary_particle_file << "\" ...\n\n";
            }
	  //                                                                                                                                                                                                     
	  // The second argument is how many Reals we read into `m_data[]`                                                                                                                                       
	  // after reading in `m_pos[]` in each of the binary particle files.                                                                                                                                    
	  // Here we're reading in the particle mass and velocity. This might 
	  // have to be extended to all Gauss Beam arguments.                                                                        
	  FDMwkbPC->InitFromBinaryMetaFile(binary_particle_file, BL_SPACEDIM + 1);
        }
      else if (particle_init_type == "BinaryMortonFile")
        {
          if (verbose)
            {
	      amrex::Print() << "\nInitializing FDM particles from morton-ordered binary file\""
                             << binary_particle_file << "\" ...\n\n";
            }
	  //                                                                                                                                                                                                     
	  // The second argument is how many Reals we read into `m_data[]`                                                                                                                                       
	  // after reading in `m_pos[]` in each of the binary particle files.                                                                                                                                    
	  // Here we're reading in the particle mass and velocity. This might                                                                                               
	  // have to be extended to all Gauss Beam arguments.                                                                                                                                          
          FDMwkbPC->InitFromBinaryMortonFile(binary_particle_file,
                                         BL_SPACEDIM + 1,
                                         particle_skip_factor);
        }
      else if (particle_init_type == "Cosmological")
        {
	  // moved to Nyx_initcosmo.cpp                                                                                                                                                                      
        }
      //do init                                                                                                                                                                                                  
      // else if(particle_init_type == "Bosonstar")
      //   {
      // 	  const int initDataDim = 4; //FIXME: depends on setup... I guess                                                                                                                                        
      // 	  const Real * dx = geom.CellSize();
      // 	  const Real cur_time = state[State_Type].curTime();
      // 	  MultiFab InitData(grids,dmap,initDataDim,0);

      // 	  for (MFIter mfi(InitData); mfi.isValid(); ++mfi)
      //       {
      // 	      RealBox gridloc = RealBox(grids[mfi.index()], geom.CellSize(), geom.ProbLo());
      // 	      const Box& box = mfi.validbox();
      // 	      const int* lo = box.loVect();
      // 	      const int* hi = box.hiVect();

      // 	      // BL_FORT_PROC_CALL(INITPART, initpart)
      // 	      initpart(&level, &cur_time,
      // 		 lo, hi,
      // 		 &initDataDim, BL_TO_FORTRAN(InitData[mfi]),
      // 		 dx, gridloc.lo(), gridloc.hi());
      //       }

 
      // 	  BoxArray baWhereNot;
      // 	  if (level < parent->initialBaLevels())
      // 	    baWhereNot = parent->initialBa(level+1);
      // 	  BoxArray myBaWhereNot(baWhereNot.size());
      // 	  for (int i=0; i < baWhereNot.size(); i++)
      // 	    myBaWhereNot.set(i, baWhereNot[i]);
      // 	  if (level < parent->initialBaLevels())
      // 	    myBaWhereNot.coarsen(parent->refRatio(level));

      // 	  if(num_particle_fdm > 0)
      // 	    FDMwkbPC->InitVarCount(InitData, num_particle_fdm, myBaWhereNot, level, parent->initialBaLevels()+1);
      // 	  else
      // 	    amrex::Error("\nNeed num_particle_fdm > 0 for InitVarCount!\n\n");

      //   }
      else if(particle_init_type == "GaussianBeams")
        {
	  if (!do_dm_particles)
	    amrex::Print() << "\n DM particles are needed for the construction of the gravitational potential!!\n\n";
	  if(num_particle_fdm > 0)
	    FDMwkbPC->InitGaussianBeams(num_particle_fdm, level, parent->initialBaLevels()+1, hbaroverm, sigma_fdm, gamma_fdm, meandens, alpha_fdm, a);
	  else
	    amrex::Error("\nNeed num_particle_fdm > 0 for InitGaussianBeams!\n\n");

	  MultiFab& Ax_new = get_level(level).get_new_data(Axion_Type);
	  Ax_new.setVal(0.);

	  // if(levelmethod[level]==GBlevel){

	  // //Define neccessary number of ghost cells                                                                                                                                              
	  // int ng = ceil(Nyx::sigma_fdm*Nyx::theta_fdm/get_level(level).Geom().CellSize()[0]);

	  // //Initialize MultiFabs                                                                                                                                                                      
	  // MultiFab fdmreal(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
	  // fdmreal.setVal(0.);
	  // MultiFab fdmimag(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
	  // fdmimag.setVal(0.);

	  // const Real cur_time  = state[State_Type].curTime();
	  // const Real a_new = get_comoving_a(cur_time);

	  // //Deposit Gaussian Beams                                                                                                                                                         
	  // if(Nyx::theFDMPC())
	  //   Nyx::theFDMPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theGhostFDMPC())
	  //   Nyx::theGhostFDMPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theVirtFDMPC())
	  //   Nyx::theVirtFDMPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theFDMwkbPC())
	  //   Nyx::theFDMwkbPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theGhostFDMwkbPC())
	  //   Nyx::theGhostFDMwkbPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theVirtFDMwkbPC())
	  //   Nyx::theVirtFDMwkbPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);

	  // //Update real part in FDM state                                                                                                                                                            
	  // Ax_new.ParallelCopy(fdmreal, 0, Nyx::AxRe, 1, fdmreal.nGrow(),
	  // 		      Ax_new.nGrow(), parent->Geom(level).periodicity(),FabArrayBase::ADD);

	  // //Update imaginary part in FDM state                                                                                                                                                   
	  // Ax_new.ParallelCopy(fdmimag, 0, Nyx::AxIm, 1, fdmimag.nGrow(),
	  // 		      Ax_new.nGrow(), parent->Geom(level).periodicity(),FabArrayBase::ADD);

	  // //Update density in FDM state                                                                                                                                                                 
	  // AmrLevel* amrlev = &parent->getLevel(level);
	  // for (amrex::FillPatchIterator fpi(*amrlev,Ax_new); fpi.isValid(); ++fpi)
	  //   {
	  //     if (Ax_new[fpi].contains_nan())
	  // 	amrex::Abort("Nans in state just before FDM density update");
	  //     BL_FORT_PROC_CALL(FORT_FDM_FIELDS, fort_fdm_fields)
	  // 	(BL_TO_FORTRAN(Ax_new[fpi]));
	  //     if (Ax_new[fpi].contains_nan())
	  // 	amrex::Abort("Nans in state just after FDM density update");
	  //   }
	  // }

        } else {
	//FIXME                                                                                                                                                                                                            
	std::cout << "multilevel init not implemented, yet!" << std::endl;
      }

    }else{
      BL_ASSERT (FDMPC == 0);
      FDMPC = new FDMParticleContainer(parent);

      if (parent->subCycle())
        {
      	  VirtFDMPC = new FDMParticleContainer(parent);

      	  GhostFDMPC = new FDMParticleContainer(parent);
        }

      //                                                                                                                                                                                                         
      // Make sure to call RemoveParticlesOnExit() on exit.                                                                                                                                                      
      //   (if do_dm_particles then we have already called ExecOnFinalize)                                                                                                                                       
      //                                                                                                                                                                                                         
      if (!do_dm_particles)
	amrex::ExecOnFinalize(RemoveParticlesOnExit);
      //                                                                                                                                                                                                         
      // 2 gives more stuff than 1.                                                                                                                                                                              
      //                                                                                                                                                                                                         
      FDMPC->SetVerbose(particle_verbose);

      if (particle_init_type == "AsciiFile")
        {
	  if (verbose && ParallelDescriptor::IOProcessor())
            {
	      amrex::Print() << "\nInitializing FDM particles from \""
			     << ascii_particle_file << "\" ...\n\n";
            }
	  //                                                                                                                                                                                                     
	  // The second argument is how many Reals we read into `m_data[]`                                                                                                                                       
	  // after reading in `m_pos[]`. Here we're reading in the particle                                                                                                                                      
	  // mass and velocity. This might have to be extended to all Gauss 
	  // Beam arguments.                                                                                                                                                               
	  FDMPC->InitFromAsciiFile(ascii_particle_file, BL_SPACEDIM + 1, &Nrep);
        }
      else if (particle_init_type == "BinaryFile")
        {
	  if (verbose)
            {
	      amrex::Print() << "\nInitializing FDM particles from \""
			     << binary_particle_file << "\" ...\n\n";
            }
	  //                                                                                                                                                                                                     
	  // The second argument is how many Reals we read into `m_data[]`                                                                                                                                       
	  // after reading in `m_pos[]`. Here we're reading in the particle                                                                                                                                      
	  // mass and velocity. This might have to be extended to all Gauss                                                                                                                          
          // Beam arguments.                                                                                                                                                                           
	  FDMPC->InitFromBinaryFile(binary_particle_file, BL_SPACEDIM + 1);
        }
      else if (particle_init_type == "BinaryMetaFile")
        {
	  if (verbose)
            {
	      amrex::Print() << "\nInitializing FDM particles from meta file\""
			     << binary_particle_file << "\" ...\n\n";
            }
	  //                                                                                                                                                                                                     
	  // The second argument is how many Reals we read into `m_data[]`                                                                                                                                       
	  // after reading in `m_pos[]` in each of the binary particle files.                                                                                                                                    
	  // Here we're reading in the particle mass and velocity. This might 
	  // have to be extended to all Gauss Beam arguments.                                                                        
	  FDMPC->InitFromBinaryMetaFile(binary_particle_file, BL_SPACEDIM + 1);
        }
      else if (particle_init_type == "BinaryMortonFile")
        {
          if (verbose)
            {
	      amrex::Print() << "\nInitializing FDM particles from morton-ordered binary file\""
                             << binary_particle_file << "\" ...\n\n";
            }
	  //                                                                                                                                                                                                     
	  // The second argument is how many Reals we read into `m_data[]`                                                                                                                                       
	  // after reading in `m_pos[]` in each of the binary particle files.                                                                                                                                    
	  // Here we're reading in the particle mass and velocity. This might                                                                                               
	  // have to be extended to all Gauss Beam arguments.                                                                                                                                          
          FDMPC->InitFromBinaryMortonFile(binary_particle_file,
                                         BL_SPACEDIM + 1,
                                         particle_skip_factor);
        }
      else if (particle_init_type == "Cosmological")
        {
	  // moved to Nyx_initcosmo.cpp                                                                                                                                                                      
        }
      //do init                                                                                                                                                                                                  
      // else if(particle_init_type == "Bosonstar")
      //   {
      // 	  const int initDataDim = 4; //FIXME: depends on setup... I guess                                                                                                                                        
      // 	  const Real * dx = geom.CellSize();
      // 	  const Real cur_time = state[State_Type].curTime();
      // 	  MultiFab InitData(grids,dmap,initDataDim,0);

      // 	  for (MFIter mfi(InitData); mfi.isValid(); ++mfi)
      //       {
      // 	      RealBox gridloc = RealBox(grids[mfi.index()], geom.CellSize(), geom.ProbLo());
      // 	      const Box& box = mfi.validbox();
      // 	      const int* lo = box.loVect();
      // 	      const int* hi = box.hiVect();

      // 	      // BL_FORT_PROC_CALL(INITPART, initpart)
      // 	      initpart(&level, &cur_time,
      // 		 lo, hi,
      // 		 &initDataDim, BL_TO_FORTRAN(InitData[mfi]),
      // 		 dx, gridloc.lo(), gridloc.hi());
      //       }

 
      // 	  BoxArray baWhereNot;
      // 	  if (level < parent->initialBaLevels())
      // 	    baWhereNot = parent->initialBa(level+1);
      // 	  BoxArray myBaWhereNot(baWhereNot.size());
      // 	  for (int i=0; i < baWhereNot.size(); i++)
      // 	    myBaWhereNot.set(i, baWhereNot[i]);
      // 	  if (level < parent->initialBaLevels())
      // 	    myBaWhereNot.coarsen(parent->refRatio(level));

      // 	  if(num_particle_fdm > 0)
      // 	    FDMPC->InitVarCount(InitData, num_particle_fdm, myBaWhereNot, level, parent->initialBaLevels()+1);
      // 	  else
      // 	    amrex::Error("\nNeed num_particle_fdm > 0 for InitVarCount!\n\n");

      //   }
      else if(particle_init_type == "GaussianBeams")
        {
	  if (!do_dm_particles)
	    amrex::Print() << "\n DM particles are needed for the construction of the gravitational potential!!\n\n";
	  if(num_particle_fdm > 0)
	    FDMPC->InitGaussianBeams(num_particle_fdm, level, parent->initialBaLevels()+1, hbaroverm, sigma_fdm, gamma_fdm, meandens, alpha_fdm, a);
	  else
	    amrex::Error("\nNeed num_particle_fdm > 0 for InitGaussianBeams!\n\n");

	  MultiFab& Ax_new = get_level(level).get_new_data(Axion_Type);
	  Ax_new.setVal(0.);
	  
	  // if(levelmethod[level]==GBlevel){
	  
	  // //Define neccessary number of ghost cells                                                                                                                                              
	  // int ng = ceil(Nyx::sigma_fdm*Nyx::theta_fdm/get_level(level).Geom().CellSize()[0]);

	  // //Initialize MultiFabs                                                                                                                                                                      
	  // MultiFab fdmreal(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
	  // fdmreal.setVal(0.);
	  // MultiFab fdmimag(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, ng);
	  // fdmimag.setVal(0.);

          // const Real cur_time  = state[State_Type].curTime();
	  // const Real a_new = get_comoving_a(cur_time);

	  // //Deposit Gaussian Beams                                                                                                                                                         
	  // if(Nyx::theFDMPC())
	  //   Nyx::theFDMPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theGhostFDMPC())
	  //   Nyx::theGhostFDMPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theVirtFDMPC())
	  //   Nyx::theVirtFDMPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theFDMwkbPC())
	  //   Nyx::theFDMwkbPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theGhostFDMwkbPC())
	  //   Nyx::theGhostFDMwkbPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);
	  // if(Nyx::theVirtFDMwkbPC())
	  //   Nyx::theVirtFDMwkbPC()->DepositFDMParticles(fdmreal,fdmimag,level,a_new);

	  // //Update real part in FDM state                                                                                                                                                            
	  // Ax_new.ParallelCopy(fdmreal, 0, Nyx::AxRe, 1, fdmreal.nGrow(),
	  // 		      Ax_new.nGrow(), parent->Geom(level).periodicity(),FabArrayBase::ADD);

	  // //Update imaginary part in FDM state                                                                                                                                                   
	  // Ax_new.ParallelCopy(fdmimag, 0, Nyx::AxIm, 1, fdmimag.nGrow(),
	  // 		      Ax_new.nGrow(), parent->Geom(level).periodicity(),FabArrayBase::ADD);

	  // //Update density in FDM state                                                                                                                                                                 
	  // AmrLevel* amrlev = &parent->getLevel(level);
	  // for (amrex::FillPatchIterator fpi(*amrlev,Ax_new); fpi.isValid(); ++fpi)
	  //   {
	  //     if (Ax_new[fpi].contains_nan())
	  // 	amrex::Abort("Nans in state just before FDM density update");
	  //     BL_FORT_PROC_CALL(FORT_FDM_FIELDS, fort_fdm_fields)
	  // 	(BL_TO_FORTRAN(Ax_new[fpi]));
	  //     if (Ax_new[fpi].contains_nan())
	  // 	amrex::Abort("Nans in state just after FDM density update");
	  //   }
	  // }
	  
        } else {
	//FIXME                                                                                                                                                                                                            
	std::cout << "multilevel init not implemented, yet!" << std::endl;
      }

    }
    }
#endif
}

#ifdef GRAVITY
#ifndef NO_HYDRO
void
Nyx::init_santa_barbara (int init_sb_vels)
{
    BL_PROFILE("Nyx::init_santa_barbara()");
    Real cur_time = state[State_Type].curTime();
    Real a = old_a;

    amrex::Print() << "... time and comoving a when data is initialized at level " 
                   << level << " " << cur_time << " " << a << '\n';

    if (level == 0)
    {
        Real frac_for_hydro = comoving_OmB / comoving_OmM;
        Real omfrac = 1.0 - frac_for_hydro;
 
        if ( (init_with_sph_particles == 0) && (frac_for_hydro != 1.0) ) {
            DMPC->MultiplyParticleMass(level, omfrac);
	}

        Vector<std::unique_ptr<MultiFab> > particle_mf(1);
        if (init_sb_vels == 1)
        {
            if (init_with_sph_particles == 1) {
               SPHPC->AssignDensityAndVels(particle_mf);
	    } else {
               DMPC->AssignDensityAndVels(particle_mf);
	    }

        } else {
            if (init_with_sph_particles == 1) {
                SPHPC->AssignDensity(particle_mf);
	    } else {
                DMPC->AssignDensity(particle_mf);
	    }
        }

        // As soon as we have used the SPH particles to define the density
        //    and velocity on the grid, we can go ahead and destroy them.
        if (init_with_sph_particles == 1) {
            delete SPHPC;
	}

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev],
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1,
                                 parent->refRatio(lev));
        }

        // Only multiply the density, not the velocities
        if (init_with_sph_particles == 0)
        {
            if (frac_for_hydro == 1.0)
            {
                particle_mf[level]->mult(0,0,1);
            }
            else
            {
                particle_mf[level]->mult(frac_for_hydro / omfrac,0,1);
            }
        }

        const Real * dx = geom.CellSize();
        MultiFab& S_new = get_new_data(State_Type);
        MultiFab& D_new = get_new_data(DiagEOS_Type);
#ifdef FDM
        MultiFab&   Ax_new   = get_new_data(Axion_Type);
        int         na       = Ax_new.nComp();
#endif

        int ns = S_new.nComp();
        int nd = D_new.nComp(); 
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            RealBox gridloc = RealBox(grids[mfi.index()],
                                      geom.CellSize(),
                                      geom.ProbLo());
            const Box& box = mfi.validbox();
            const int* lo = box.loVect();
            const int* hi = box.hiVect();

            // Temp unused for GammaLaw, set it here so that pltfiles have
            // defined numbers
            D_new[mfi].setVal(0, Temp_comp);
            D_new[mfi].setVal(0,   Ne_comp);

            fort_initdata
                (level, cur_time, lo, hi, 
                 ns,BL_TO_FORTRAN(S_new[mfi]), 
#ifdef FDM
                 na, BL_TO_FORTRAN(Ax_new[mfi]),
#endif
                 nd,BL_TO_FORTRAN(D_new[mfi]), dx,
                 gridloc.lo(), gridloc.hi());
        }

        if (inhomo_reion) init_zhi();

        // Add the particle density to the gas density 
        MultiFab::Add(S_new, *particle_mf[level], 0, Density, 1, S_new.nGrow());

        if (init_sb_vels == 1)
        {
            // Convert velocity to momentum
            for (int i = 0; i < BL_SPACEDIM; ++i) {
               MultiFab::Multiply(*particle_mf[level], *particle_mf[level], 0, 1+i, 1, 0);
	    }

            // Add the particle momenta to the gas momenta (initially zero)
            MultiFab::Add(S_new, *particle_mf[level], 1, Xmom, BL_SPACEDIM, S_new.nGrow());
        }

    } else {

        MultiFab& S_new = get_new_data(State_Type);
        FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, S_new.nComp());

        MultiFab& D_new = get_new_data(DiagEOS_Type);
        FillCoarsePatch(D_new, 0, cur_time, DiagEOS_Type, 0, D_new.nComp());
 
        MultiFab& Phi_new = get_new_data(PhiGrav_Type);
        FillCoarsePatch(Phi_new, 0, cur_time, PhiGrav_Type, 0, Phi_new.nComp());

       // Convert (rho X)_i to X_i before calling init_e_from_T
       if (use_const_species == 0) {
           for (int i = 0; i < NumSpec; ++i) {
               MultiFab::Divide(S_new, S_new, Density, FirstSpec+i, 1, 0);
	   }
       }
    }

    // Make sure we've finished initializing the density before calling this.
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);
    int ns = S_new.nComp();
    int nd = D_new.nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        const int* lo = box.loVect();
        const int* hi = box.hiVect();

        fort_init_e_from_t
            (BL_TO_FORTRAN(S_new[mfi]), &ns, 
             BL_TO_FORTRAN(D_new[mfi]), &nd, lo, hi, &a);
    }

    // Convert X_i to (rho X)_i
    if (use_const_species == 0) {
        for (int i = 0; i < NumSpec; ++i) {
            MultiFab::Multiply(S_new, S_new, Density, FirstSpec+i, 1, 0);
	}
    }
}
#endif
#endif

void
Nyx::particle_post_restart (const std::string& restart_file, bool is_checkpoint)
{
    BL_PROFILE("Nyx::particle_post_restart()");

    if (level > 0)
        return;
     
    if (do_dm_particles)
    {
        BL_ASSERT(DMPC == 0);
        DMPC = new DarkMatterParticleContainer(parent);
        ActiveParticles.push_back(DMPC);
 
        if (parent->subCycle())
        {
            VirtPC = new DarkMatterParticleContainer(parent);
            VirtualParticles.push_back(VirtPC);
 
            GhostPC = new DarkMatterParticleContainer(parent);
            GhostParticles.push_back(GhostPC);
        }

        //
        // Make sure to call RemoveParticlesOnExit() on exit.
        //
        amrex::ExecOnFinalize(RemoveParticlesOnExit);
        //
        // 2 gives more stuff than 1.
        //
        DMPC->SetVerbose(particle_verbose);
        DMPC->Restart(restart_file, dm_chk_particle_file, is_checkpoint);
        //
        // We want the ability to write the particles out to an ascii file.
        //
        ParmParse pp("particles");

        std::string dm_particle_output_file;
        pp.query("dm_particle_output_file", dm_particle_output_file);

        if (!dm_particle_output_file.empty())
        {
            DMPC->WriteAsciiFile(dm_particle_output_file);
        }
    }
#ifdef AGN
    {
        BL_ASSERT(APC == 0);
        APC = new AGNParticleContainer(parent, num_particle_ghosts);
        ActiveParticles.push_back(APC);
 
        if (parent->subCycle())
        {
          VirtAPC = new AGNParticleContainer(parent, num_particle_ghosts);
          VirtualParticles.push_back(VirtAPC);
 
          GhostAPC = new AGNParticleContainer(parent, num_particle_ghosts);
          GhostParticles.push_back(GhostAPC);
        }

        //
        // Make sure to call RemoveParticlesOnExit() on exit.
        //
        amrex::ExecOnFinalize(RemoveParticlesOnExit);
        //
        // 2 gives more stuff than 1.
        //
        APC->SetVerbose(particle_verbose);
        APC->Restart(restart_file, agn_chk_particle_file, is_checkpoint);
        //
        // We want the ability to write the particles out to an ascii file.
        //
        ParmParse pp("particles");

        std::string agn_particle_output_file;
        pp.query("agn_particle_output_file", agn_particle_output_file);

        if (!agn_particle_output_file.empty())
        {
            APC->WriteAsciiFile(agn_particle_output_file);
        }
    }
#endif
#ifdef FDM
    {
      if(partlevel){
      if(wkb_approx){
	BL_ASSERT (FDMwkbPC == 0);
	FDMwkbPC = new FDMwkbParticleContainer(parent);
	
	// ActiveParticles.push_back(FDMwkbPC);
	
	if (parent->subCycle())
	  {
	    VirtFDMwkbPC = new FDMwkbParticleContainer(parent);
	    // VirtualParticles.push_back(VirtFDMwkbPC);
	    
	    GhostFDMwkbPC = new FDMwkbParticleContainer(parent);
	    // GhostParticles.push_back(GhostFDMwkbPC);
	  }
	
	//                                                                                                                                                                                                     
	// Make sure to call RemoveParticlesOnExit() on exit.                                                                                                                                                  
	//   (if do_dm_particles then we have already called ExecOnFinalize)                                                                                                                                  
	//                                                                                                                                                                                                   
	if (!do_dm_particles)
	  amrex::ExecOnFinalize(RemoveParticlesOnExit);
	//                                                                                                                                                                                                      
	// 2 gives more stuff than 1.                                                                                                                                                                          
	//                                                                                                                                                                                                     
	FDMwkbPC->SetVerbose(particle_verbose);
	FDMwkbPC->Restart(restart_file, fdm_chk_particle_file, is_checkpoint);
	//                                                                                                                                                                                                       
	// We want the ability to write the particles out to an ascii file.                                                                                                                                     
	//                                                                                                                                                                                                      
	ParmParse pp("particles");
	
	std::string fdm_particle_output_file;
	pp.query("fdm_particle_output_file", fdm_particle_output_file);
	
	if (!fdm_particle_output_file.empty())
	  {
	    FDMwkbPC->WriteAsciiFile(fdm_particle_output_file);
	  }
      }else{
	BL_ASSERT (FDMPC == 0);
	FDMPC = new FDMParticleContainer(parent);
	
	// ActiveParticles.push_back(FDMPC);
	
	if (parent->subCycle())
	  {
	    VirtFDMPC = new FDMParticleContainer(parent);
	    // VirtualParticles.push_back(VirtFDMPC);
	    
	    GhostFDMPC = new FDMParticleContainer(parent);
	    // GhostParticles.push_back(GhostFDMPC);
	  }
	
	//                                                                                                                                                                                                     
	// Make sure to call RemoveParticlesOnExit() on exit.                                                                                                                                                  
	//   (if do_dm_particles then we have already called ExecOnFinalize)                                                                                                                                   
	//                                                                                                                                                                                                     
	if (!do_dm_particles)
	  amrex::ExecOnFinalize(RemoveParticlesOnExit);
	//                                                                                                                                                                                                     
	// 2 gives more stuff than 1.                                                                                                                                                                          
	//                                                                                                                                                                                                     
	FDMPC->SetVerbose(particle_verbose);
	FDMPC->Restart(restart_file, fdm_chk_particle_file, is_checkpoint);
	//                                                                                                                                                                                                    
	// We want the ability to write the particles out to an ascii file.                                                                                                                                   
	//                                                                                                                                                                                                    
	ParmParse pp("particles");
	
	std::string fdm_particle_output_file;
	pp.query("fdm_particle_output_file", fdm_particle_output_file);
	
	if (!fdm_particle_output_file.empty())
	  {
	    FDMPC->WriteAsciiFile(fdm_particle_output_file);
	  }
      }
    }
    }
#endif
}

#ifdef GRAVITY
void
Nyx::particle_est_time_step (Real& est_dt)
{
    BL_PROFILE("Nyx::particle_est_time_step()");

    const Real cur_time = state[PhiGrav_Type].curTime();
    const Real a = get_comoving_a(cur_time);
    MultiFab& grav = get_new_data(Gravity_Type);

    if (DMPC && particle_move_type == "Gravitational")
      {
        const Real est_dt_particle = DMPC->estTimestep(grav, a, level, particle_cfl);
	
        if (est_dt_particle > 0) {
	  est_dt = std::min(est_dt, est_dt_particle);
	}

        if (verbose)
        {
            if (est_dt_particle > 0)
            {
                amrex::Print() << "...estdt from particles at level "
                          << level << ": " << est_dt_particle << '\n';
            }
            else
            {
                amrex::Print() << "...there are no particles at level "
                          << level << '\n';
            }
	}
      }
#ifdef NEUTRINO_PARTICLES
    if (NPC && particle_move_type == "Gravitational")
      {
        const Real est_dt_neutrino = NPC->estTimestep(grav, a, level, neutrino_cfl);

        if (est_dt_neutrino > 0) {
	  est_dt = std::min(est_dt, est_dt_neutrino);
	}

        if (verbose)
        {
            if (est_dt_neutrino > 0)
            {
                amrex::Print() << "...estdt from neutrinos at level "
                          << level << ": " << est_dt_neutrino << '\n';
            }
            else
            {
                amrex::Print() << "...there are no neutrinos at level "
                          << level << '\n';
            }
	}
      }
#endif

#ifdef FDM
    if (FDMwkbPC && particle_move_type == "Gravitational")
      {
	// const Real est_dt_fdm1 = FDMwkbPC->estTimestep(grav, a, level, particle_cfl);
	// MultiFab& phi = get_new_data(PhiGrav_Type);
	// Real beam_cfl_reduced = beam_cfl*2.0*M_PI*hbaroverm;
	// const Real est_dt_fdm2 = FDMwkbPC->estTimestepFDM(phi, level, beam_cfl_reduced);
	// const Real est_dt_fdm = std::min(est_dt_fdm1 ,est_dt_fdm2);
	const Real est_dt_fdm = FDMwkbPC->estTimestep(grav, a, level, particle_cfl);
	
	if (est_dt_fdm > 0)
	  est_dt = std::min(est_dt, est_dt_fdm);
	
	if (verbose)
	  {
	    if (est_dt_fdm > 0)
	      {
		amrex::Print() << "...estdt from FDMwkb at level "
			       << level << ": " << est_dt_fdm << '\n';
	      }
	    else
	      {
		amrex::Print() << "...there are no FDMwkb beams at level "
			       << level << '\n';
	      }
	  }
      }
    if (FDMPC && particle_move_type == "Gravitational")
      {
	const Real est_dt_fdm1 = FDMPC->estTimestep(grav, a, level, particle_cfl);
	MultiFab& phi = get_new_data(PhiGrav_Type);
	Real beam_cfl_reduced = beam_cfl*2.0*M_PI*hbaroverm;
	const Real est_dt_fdm2 = FDMPC->estTimestepFDM(phi, level, beam_cfl_reduced);
	const Real est_dt_fdm = std::min(est_dt_fdm1 ,est_dt_fdm2);
	// const Real est_dt_fdm = FDMPC->estTimestep(grav, a, level, particle_cfl);
	
	if (est_dt_fdm > 0)
	  est_dt = std::min(est_dt, est_dt_fdm);
	  
	if (verbose)
	  {
	    if (est_dt_fdm > 0)
	      {
		amrex::Print() << "...estdt from FDM at level "
			       << level << ": " << est_dt_fdm << '\n';
	      }
	    else
	      {
		amrex::Print() << "...there are no FDM beams at level "
			       << level << '\n';
	      }
	  }
      }
#endif
}
#endif

void
Nyx::particle_redistribute (int lbase, bool my_init)
{
    BL_PROFILE("Nyx::particle_redistribute()");
#ifdef FDM
    if (DMPC || FDMPC || FDMwkbPC)
#else
    if (DMPC)
#endif
    {
        //  
        // If we are calling with my_init = true, then we want to force the redistribute
        //    without checking whether the grids have changed.
        //  
        if (my_init)
        {
#ifdef FDM
	    if(DMPC)
	      DMPC->Redistribute(lbase);
	    if(FDMPC)
	      FDMPC->Redistribute(lbase);
	    if(FDMwkbPC)
	      FDMwkbPC->Redistribute(lbase);
#else
            DMPC->Redistribute(lbase);
#endif
            return;
        }

        //
        // These are usually the BoxArray and DMap from the last regridding.
        //
        static Vector<BoxArray>            ba;
        static Vector<DistributionMapping> dm;

        bool    changed      = false;
        bool dm_changed      = false;
        bool ba_changed      = false;
        bool ba_size_changed = false;

        int flev = parent->finestLevel();
	
        while ( parent->getAmrLevels()[flev] == nullptr ) {
            flev--;
	}
 
        if (ba.size() != flev+1)
        {
            amrex::Print() << "BA SIZE " << ba.size() << std::endl;
            amrex::Print() << "FLEV    " << flev << std::endl;
            ba.resize(flev+1);
            dm.resize(flev+1);
            changed = true;
            ba_size_changed = true;
        }
        else
        {
            for (int i = 0; i <= flev && !changed; i++)
            {
                if (ba[i] != parent->boxArray(i))
                {
                    //
                    // The BoxArrays have changed in the regridding.
                    //
                    changed = true;
                    ba_changed = true;
                }

                if ( ! changed)
                {
                    if (dm[i] != parent->getLevel(i).get_new_data(0).DistributionMap())
                    {
                        //
                        // The DistributionMaps have changed in the regridding.
                        //
                        changed = true;
                        dm_changed = true;
                    }
                }
            }
        }

        if (changed)
        {
            //
            // We only need to call Redistribute if the BoxArrays or DistMaps have changed.
	    // We also only call it for particles >= lbase. This is
	    // because of we called redistribute during a subcycle, there may be particles not in
	    // the proper position on coarser levels.
            //
            if (verbose)
            {
                if (ba_size_changed)
                   amrex::Print() << "Calling redistribute because the size of BoxArray changed " << '\n';
                else if (ba_changed)
                   amrex::Print() << "Calling redistribute because BoxArray changed " << '\n';
                else if (dm_changed)
                   amrex::Print() << "Calling redistribute because DistMap changed " << '\n';
            }

#ifdef FDM
	    if(DMPC)
	      DMPC->Redistribute(lbase);
	    if(FDMPC)
	      FDMPC->Redistribute(lbase);
	    if(FDMwkbPC)
	      FDMwkbPC->Redistribute(lbase);
#else
            DMPC->Redistribute(lbase);
#endif
            //
            // Use the new BoxArray and DistMap to define ba and dm for next time.
            //
            for (int i = 0; i <= flev; i++)
            {
                ba[i] = parent->boxArray(i);
                dm[i] = parent->getLevel(i).get_new_data(0).DistributionMap();
            }
        }
        else
        {
            if (verbose)
                amrex::Print() << "NOT calling redistribute because NOT changed " << '\n';
        }
    }
}

void
Nyx::particle_move_random ()
{
    BL_PROFILE("Nyx::particle_move_random()");
    if (DMPC && particle_move_type == "Random")
    {
        BL_ASSERT(level == 0);

        DMPC->MoveRandom();
    }
}

void
Nyx::setup_virtual_particles()
{
  BL_PROFILE("Nyx::setup_virtual_particles()");
  if(!virtual_particles_set)
    {
      if(Nyx::theDMPC() != 0){
        DarkMatterParticleContainer::AoS virts;
        if (level < parent->finestLevel())
	  {
    	    get_level(level + 1).setup_virtual_particles();
	    Nyx::theVirtPC()->CreateVirtualParticles(level+1, virts);
	    Nyx::theVirtPC()->AddParticlesAtLevel(virts, level);
	    Nyx::theDMPC()->CreateVirtualParticles(level+1, virts);
	    Nyx::theVirtPC()->AddParticlesAtLevel(virts, level);
	  }
      }
#ifdef FDM
      if(Nyx::theFDMPC() != 0){
        FDMParticleContainer::AoS virts;
        if (level < parent->finestLevel())
	  {
	    if(Nyx::theDMPC() == 0)
	      get_level(level + 1).setup_virtual_particles();
	    Nyx::theVirtFDMPC()->CreateVirtualParticles(level+1, virts);
	    Nyx::theVirtFDMPC()->AddParticlesAtLevel(virts, level);
	    Nyx::theFDMPC()->CreateVirtualParticles(level+1, virts);
	    Nyx::theVirtFDMPC()->AddParticlesAtLevel(virts, level);
	  }
      }
      if(Nyx::theFDMwkbPC() != 0){
        FDMwkbParticleContainer::AoS virts;
        if (level < parent->finestLevel())
	  {
	    if(Nyx::theDMPC() == 0 && Nyx::theFDMPC() == 0)
	      get_level(level + 1).setup_virtual_particles();
	    Nyx::theVirtFDMwkbPC()->CreateVirtualParticles(level+1, virts);
	    Nyx::theVirtFDMwkbPC()->AddParticlesAtLevel(virts, level);
	    Nyx::theFDMwkbPC()->CreateVirtualParticles(level+1, virts);
	    Nyx::theVirtFDMwkbPC()->AddParticlesAtLevel(virts, level);
	  }
      }
#endif
      virtual_particles_set = true;
    }
}

void
Nyx::remove_virtual_particles()
{
    BL_PROFILE("Nyx::remove_virtual_particles()");
    for (int i = 0; i < VirtualParticles.size(); i++)
    {
        if (VirtualParticles[i] != 0)
            VirtualParticles[i]->RemoveParticlesAtLevel(level);
        virtual_particles_set = false;
    }
#ifdef FDM
    if(theVirtFDMPC()){
      theVirtFDMPC()->RemoveParticlesAtLevel(level);
      virtual_particles_set = false;
    }
    if(theVirtFDMwkbPC()){
      theVirtFDMwkbPC()->RemoveParticlesAtLevel(level);
      virtual_particles_set = false;
    }
#endif
}

void
Nyx::setup_ghost_particles(int ngrow)
{
    BL_PROFILE("Nyx::setup_ghost_particles()");
    BL_ASSERT(level < parent->finestLevel());
    if(Nyx::theDMPC() != 0)
    {
        DarkMatterParticleContainer::AoS ghosts;
        Nyx::theDMPC()->CreateGhostParticles(level, ngrow, ghosts);
        Nyx::theGhostPC()->AddParticlesAtLevel(ghosts, level+1, ngrow);
    }
#ifdef AGN
    if(Nyx::theAPC() != 0)
    {
        AGNParticleContainer::AoS ghosts;
        Nyx::theAPC()->CreateGhostParticles(level, ngrow, ghosts);
        Nyx::theGhostAPC()->AddParticlesAtLevel(ghosts, level+1, ngrow);
    }
#endif
#ifdef NEUTRINO_PARTICLES
    if(Nyx::theNPC() != 0)
    {
        NeutrinoParticleContainer::AoS ghosts;
        Nyx::theNPC()->CreateGhostParticles(level, ngrow, ghosts);
        Nyx::theGhostNPC()->AddParticlesAtLevel(ghosts, level+1, ngrow);
    }
#endif
#ifdef FDM
    // if(Nyx::theFDMPC() != 0)
    // {
    //     FDMParticleContainer::AoS ghosts;
    // 	int ng = ceil(sigma_fdm*theta_fdm/parent->Geom(level+1).CellSize()[0]);
    //     Nyx::theFDMPC()->CreateGhostParticles(level, ng, ghosts);
    //     Nyx::theGhostFDMPC()->AddParticlesAtLevel(ghosts, level+1, ng);
    // }
    // if(Nyx::theFDMwkbPC() != 0)
    // {
    //     FDMwkbParticleContainer::AoS ghosts;
    //     Nyx::theFDMwkbPC()->CreateGhostParticles(level, ngrow, ghosts);
    //     Nyx::theGhostFDMwkbPC()->AddParticlesAtLevel(ghosts, level+1, ngrow);
    // }
    if(Nyx::theFDMPC() != 0)
      {
    	FDMParticleContainer::AoS ghosts;
    	for (int lev = level+1; lev <= parent->finestLevel(); lev++){
    	  if(levelmethod[lev]==GBlevel){
    	    int ng = parent->nCycle(level)+ceil(sigma_fdm*theta_fdm/parent->Geom(lev).CellSize()[0]);
    	    Nyx::theFDMPC()->CreateGhostParticlesFDM(level, lev, ng, ghosts);
    	    Nyx::theGhostFDMPC()->AddParticlesAtLevel(ghosts, lev, ng);
    	  }
    	}
      }
    if(Nyx::theFDMwkbPC() != 0)
      {
    	FDMwkbParticleContainer::AoS ghosts;
    	for (int lev = level+1; lev <= parent->finestLevel(); lev++){
    	  if(levelmethod[lev]==GBlevel){
    	    int ng = parent->nCycle(level)+ceil(sigma_fdm*theta_fdm/parent->Geom(lev).CellSize()[0]);
    	    Nyx::theFDMwkbPC()->CreateGhostParticlesFDM(level, lev, ng, ghosts);
    	    Nyx::theGhostFDMwkbPC()->AddParticlesAtLevel(ghosts, lev, ng);
    	  }
    	}
      }
#endif

}

void
Nyx::remove_ghost_particles()
{
    BL_PROFILE("Nyx::remove_ghost_particles()");
    for (int i = 0; i < GhostParticles.size(); i++)
    {
        if (GhostParticles[i] != 0)
            GhostParticles[i]->RemoveParticlesAtLevel(level);
    }
#ifdef FDM
    if(theGhostFDMPC())
      theGhostFDMPC()->RemoveParticlesAtLevel(level);
    if(theGhostFDMwkbPC())
      theGhostFDMwkbPC()->RemoveParticlesAtLevel(level);
#endif
}

//NyxParticleContainerBase::~NyxParticleContainerBase() {}
