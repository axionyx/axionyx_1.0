#ifdef FDM
#include <stdint.h>

#include "FDMParticleContainer.H"
#include "fdm_F.H"

using namespace amrex;

/// These are helper functions used when initializing from a morton-ordered
/// binary particle file.
namespace {

  inline uint64_t split(unsigned int a) {
    uint64_t x = a & 0x1fffff;
    x = (x | x << 32) & 0x1f00000000ffff;
    x = (x | x << 16) & 0x1f0000ff0000ff;
    x = (x | x << 8)  & 0x100f00f00f00f00f;
    x = (x | x << 4)  & 0x10c30c30c30c30c3;
    x = (x | x << 2)  & 0x1249249249249249;
    return x;
  }
  
  inline uint64_t get_morton_index(unsigned int x,
				   unsigned int y,
				   unsigned int z) {
    uint64_t morton_index = 0;
    morton_index |= split(x) | ( split(y) << 1) | (split(z) << 2);
    return morton_index;
  }  

  struct BoxMortonKey {
    uint64_t morton_id;
    int box_id;
  };

  struct by_morton_id { 
    bool operator()(const BoxMortonKey &a, const BoxMortonKey &b) { 
      return a.morton_id < b.morton_id;
    }
  };

  std::string get_file_name(const std::string& base, int file_num) {
    std::stringstream ss;
    ss << base << file_num;
    return ss.str();
  }

  struct ParticleMortonFileHeader {
    long NP;
    int  DM;
    int  NX;
    int  SZ;
    int  NF;
  };
  
  void ReadHeader(const std::string& dir,
		  const std::string& file,
		  ParticleMortonFileHeader& hdr) {
    std::string header_filename = dir;
    header_filename += "/";
    header_filename += file;
    
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(header_filename, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream HdrFile(fileCharPtrString, std::istringstream::in);

    HdrFile >> hdr.NP;
    HdrFile >> hdr.DM;
    HdrFile >> hdr.NX;
    HdrFile >> hdr.SZ;
    HdrFile >> hdr.NF;    
  }

}

void
FDMParticleContainer::moveKickDrift (amrex::MultiFab&       acceleration,
		                            int                    lev,
                    			    amrex::Real            dt,
		                	    amrex::Real            a_old,
					    amrex::Real            a_half,
					    int                    where_width)
{
    BL_PROFILE("FDMParticleContainer::moveKickDrift()");

    //If there are no particles at this level
    if (lev >= this->GetParticles().size())
        return;

    const Real* dx = Geom(lev).CellSize();

    amrex::MultiFab* ac_ptr;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_ptr = &acceleration;
    }
    else
    {
        ac_ptr = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev),
			             this->m_gdb->ParticleDistributionMap(lev),
				     acceleration.nComp(),acceleration.nGrow());
        for (amrex::MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary();
    }

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 1;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();

           update_fdm_particles(&Np, particles.data(),
	   			(*ac_ptr)[pti].dataPtr(),
	   			ac_box.loVect(), ac_box.hiVect(),
	   			plo,dx,dt,a_old,a_half,&do_move);
        }
    }

    if (ac_ptr != &acceleration) delete ac_ptr;
    
    ParticleLevel&    pmap          = this->GetParticles(lev);
    if (lev > 0 && sub_cycle)
    {
        amrex::ParticleLocData pld; 
        for (auto& kv : pmap) {
            AoS&  pbox       = kv.second.GetArrayOfStructs();
            const int   n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for private(pld)
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.id() <= 0) continue;

                // Move the particle to the proper ghost cell. 
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!this->Where(p, pld, lev, lev, where_width))
                {
                    // Assert that the particle being removed is a ghost particle;
                    // the ghost particle is no longer in relevant ghost cells for this grid.
                    if (p.id() == amrex::GhostParticleID)
                    {
                        p.id() = -1;
                    }
                    else
                    {
                        std::cout << "Oops -- removing particle " << p.id() << std::endl;
                        amrex::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
        }
    }
}

void
FDMParticleContainer::moveKick (MultiFab&       acceleration,
                                       int             lev,
                                       Real            dt,
                                       Real            a_new,
                                       Real            a_half) 
{
    BL_PROFILE("FDMParticleContainer::moveKick()");

    const Real* dx = Geom(lev).CellSize();

    MultiFab* ac_ptr;
    if (OnSameGrids(lev,acceleration))
    {
        ac_ptr = &acceleration;
    }
    else 
    {
        ac_ptr = new MultiFab(ParticleBoxArray(lev),
				  ParticleDistributionMap(lev),
				  acceleration.nComp(),acceleration.nGrow());
        for (MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary();
    }

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();

           update_fdm_particles(&Np, particles.data(),
	   			(*ac_ptr)[pti].dataPtr(),
	   			ac_box.loVect(), ac_box.hiVect(),
	   			plo,dx,dt,a_half,a_new,&do_move);
        }
    }
    
    if (ac_ptr != &acceleration) delete ac_ptr;
}

void
FDMParticleContainer::moveKickDriftFDM (amrex::MultiFab&       phi,
					int                    grav_n_grow,
					amrex::MultiFab&       acceleration,
					int                    lev,
					amrex::Real            dt,
					amrex::Real            a_old,
					amrex::Real            a_half,
					int                    where_width)
{
    BL_PROFILE("FDMParticleContainer::moveKickDrift()");

    //If there are no particles at this level
    if (lev >= this->GetParticles().size())
        return;

    const Real* dx = Geom(lev).CellSize();

    amrex::MultiFab* ac_ptr;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_ptr = &acceleration;
    }
    else
    {
        ac_ptr = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev),
			             this->m_gdb->ParticleDistributionMap(lev),
				     acceleration.nComp(),acceleration.nGrow());
        for (amrex::MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary();
    }

    amrex::MultiFab* phi_ptr;
    phi_ptr = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev),
				  this->m_gdb->ParticleDistributionMap(lev),
				  phi.nComp(),grav_n_grow);
    for (amrex::MFIter mfi(*phi_ptr); mfi.isValid(); ++mfi)
      phi_ptr->setVal(0.);
    phi_ptr->copy(phi,0,0,phi.nComp());
    phi_ptr->FillBoundary();

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 1;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();
           const Box& phi_box = (*phi_ptr)[pti].box();

           update_gaussian_beams(&Np, particles.data(),
	   			 (*ac_ptr)[pti].dataPtr(),
	   			 ac_box.loVect(), ac_box.hiVect(),
	   			 (*phi_ptr)[pti].dataPtr(),
	   			 phi_box.loVect(), phi_box.hiVect(),
	   			 plo,dx,dt,a_old,a_half,&do_move);
        }
    }

    if (ac_ptr != &acceleration) delete ac_ptr;
    delete phi_ptr;
    
    ParticleLevel&    pmap          = this->GetParticles(lev);
    if (lev > 0 && sub_cycle)
    {
        amrex::ParticleLocData pld; 
        for (auto& kv : pmap) {
            AoS&  pbox       = kv.second.GetArrayOfStructs();
            const int   n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for private(pld)
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.id() <= 0) continue;

                // Move the particle to the proper ghost cell. 
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!this->Where(p, pld, lev, lev, where_width))
                {
                    // Assert that the particle being removed is a ghost particle;
                    // the ghost particle is no longer in relevant ghost cells for this grid.
                    if (p.id() == amrex::GhostParticleID)
                    {
                        p.id() = -1;
                    }
                    else
                    {
                        std::cout << "Oops -- removing particle " << p.id() << std::endl;
                        amrex::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
        }
    }
}

void
FDMParticleContainer::moveKickFDM (amrex::MultiFab& phi,
				   int              grav_n_grow,
				   amrex::MultiFab& acceleration,
				   int              lev,
				   Real             dt,
				   Real             a_new,
				   Real             a_half) 
{
    BL_PROFILE("FDMParticleContainer::moveKick()");

    const Real* dx = Geom(lev).CellSize();

    MultiFab* ac_ptr;
    if (OnSameGrids(lev,acceleration))
    {
        ac_ptr = &acceleration;
    }
    else 
    {
        ac_ptr = new MultiFab(ParticleBoxArray(lev),
				  ParticleDistributionMap(lev),
				  acceleration.nComp(),acceleration.nGrow());
        for (MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary();
    }

    MultiFab* phi_ptr;
    phi_ptr = new MultiFab(ParticleBoxArray(lev),
			   ParticleDistributionMap(lev),
			   phi.nComp(),grav_n_grow);
    for (MFIter mfi(*phi_ptr); mfi.isValid(); ++mfi)
      phi_ptr->setVal(0.);
    phi_ptr->copy(phi,0,0,phi.nComp());
    phi_ptr->FillBoundary();

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();
           const Box& phi_box = (*phi_ptr)[pti].box();

           update_gaussian_beams(&Np, particles.data(),
	   			 (*ac_ptr)[pti].dataPtr(),
	   			 ac_box.loVect(), ac_box.hiVect(),
	   			 (*phi_ptr)[pti].dataPtr(),
	   			 phi_box.loVect(), phi_box.hiVect(),
	   			 plo,dx,dt,a_half,a_new,&do_move);
        }
    }
    
    if (ac_ptr != &acceleration) delete ac_ptr;
    delete phi_ptr;
}

void
FDMParticleContainer::CreateGhostParticlesFDM (int level, int lev, int nGrow, AoS& ghosts) const
{
  BL_PROFILE("FDMParticleContainer::CreateGhostParticlesFDM()");
  BL_ASSERT(ghosts.empty());
  BL_ASSERT(level < finestLevel());

  if (lev > static_cast<int>(GetParticles().size()))
    return;

  const BoxArray& fine = ParticleBoxArray(lev);

  std::vector< std::pair<int,Box> > isects;

  const auto& pmap = GetParticles(level);
  for (const auto& kv : pmap)
    {
      const auto& pbox = kv.second.GetArrayOfStructs();
      for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
        {
	  const IntVect& iv = Index(*it, lev);
	  fine.intersections(Box(iv,iv),isects,true,nGrow);
	  if(!isects.empty())
            {
	      ParticleType p = *it;  // yes, make a copy                                                                                                                                                         
	      p.m_idata.id = GhostParticleID;
	      ghosts().push_back(p);
            }
        }
    }
}

/*
  Particle deposition
*/

void                                                                                                                                                                                                            
FDMParticleContainer::DepositFDMParticles(MultiFab& mf_real, MultiFab& mf_imag, int lev, amrex::Real a) const
{
  BL_PROFILE("FDMParticleContainer::DepositFDMParticles()");

  MultiFab* mf_pointer_real;

  if (OnSameGrids(lev, mf_real)) {
    // If we are already working with the internal mf defined on the                                                                                                                                             
    // particle_box_array, then we just work with this.                                                                                                                                                          
    mf_pointer_real = &mf_real;
  }
  else {

      amrex::Print() <<"FDM:: mf_real different!!\n";

    // If mf_real is not defined on the particle_box_array, then we need                                                                                                                                 
    // to make a temporary here and copy into mf_real at the end.                                                                                                                                        
    mf_pointer_real = new MultiFab(ParticleBoxArray(lev),
				   ParticleDistributionMap(lev),
				   1, mf_real.nGrow());

    for (MFIter mfi(*mf_pointer_real); mfi.isValid(); ++mfi) {
      (*mf_pointer_real)[mfi].setVal(0);
    }
  }

  MultiFab* mf_pointer_imag;

  if (OnSameGrids(lev, mf_imag)) {
    // If we are already working with the internal mf defined on the                                                                                                                                             
    // particle_box_array, then we just work with this.                                                                                                                                                          
    mf_pointer_imag = &mf_imag;
  }
  else {

      amrex::Print() <<"FDM:: mf_imag different!!\n";

    // If mf_imag is not defined on the particle_box_array, then we need                                                                                                                                 
    // to make a temporary here and copy into mf_imag at the end.                                                                                                                                        
      mf_pointer_imag = new MultiFab(ParticleBoxArray(lev),
				     ParticleDistributionMap(lev),
				     1, mf_imag.nGrow());

    for (MFIter mfi(*mf_pointer_imag); mfi.isValid(); ++mfi) {
      (*mf_pointer_imag)[mfi].setVal(0);
    }
  }

#ifdef _OPENMP
  const int       ng_real     = mf_pointer_real->nGrow();
  const int       ng_imag     = mf_pointer_imag->nGrow();
#endif
  const Real      strttime    = amrex::second();
  const Geometry& gm          = Geom(lev);
  const Real*     plo         = gm.ProbLo();
  const Real*     dx          = gm.CellSize();
   
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    FArrayBox local_real, local_imag;
    for (MyConstParIter pti(*this, lev); pti.isValid(); ++pti) {
      const auto& particles = pti.GetArrayOfStructs();
      const long np = pti.numParticles();
      FArrayBox& fab_real = (*mf_pointer_real)[pti];
      FArrayBox& fab_imag = (*mf_pointer_imag)[pti];
      Real *data_ptr_real, *data_ptr_imag;
      const int *lo_real, *hi_real, *lo_imag, *hi_imag;

#ifdef _OPENMP
      Box tile_box_real = pti.tilebox();
      tile_box_real.grow(ng_real);
      local_real.resize(tile_box_real,1);
      local_real = 0.0;
      data_ptr_real = local_real.dataPtr();
      lo_real = tile_box_real.loVect();
      hi_real = tile_box_real.hiVect();

      Box tile_box_imag = pti.tilebox();
      tile_box_imag.grow(ng_imag);
      local_imag.resize(tile_box_imag,1);
      local_imag = 0.0;
      data_ptr_imag = local_imag.dataPtr();
      lo_imag = tile_box_imag.loVect();
      hi_imag = tile_box_imag.hiVect();
#else
      const Box& box_real = fab_real.box();
      data_ptr_real = fab_real.dataPtr();
      lo_real = box_real.loVect();
      hi_real = box_real.hiVect();

      const Box& box_imag = fab_imag.box();
      data_ptr_imag = fab_imag.dataPtr();
      lo_imag = box_imag.loVect();
      hi_imag = box_imag.hiVect();

#endif

      //Boundaries could be easily generalized but don't need to
      BL_ASSERT(lo_real[0]==lo_imag[0]);
      BL_ASSERT(lo_real[1]==lo_imag[1]);
      BL_ASSERT(lo_real[2]==lo_imag[2]);
      BL_ASSERT(hi_real[0]==hi_imag[0]);
      BL_ASSERT(hi_real[1]==hi_imag[1]);
      BL_ASSERT(hi_real[2]==hi_imag[2]);

      deposit_fdm_particles(particles.data(), &np, data_ptr_real,
			    lo_real, hi_real, data_ptr_imag, lo_imag, 
			    hi_imag, plo, dx, a);

// #ifdef _OPENMP
//     amrex::Print() << "amrex_atomic_accumulate_fab \n";
//       amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_real),
// 				  BL_TO_FORTRAN_3D(fab_real), 1);
//       amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_imag),
// 				  BL_TO_FORTRAN_3D(fab_imag), 1);
// #endif
    }
  }

  // If mf_real is not defined on the particle_box_array, then we need                                                                                                                                   
  // to copy here from mf_pointer_real into mf_real. I believe that we don't                                                                                                                                  
  // need any information in ghost cells so we don't copy those.                                                                                                                                                 
  if (mf_pointer_real != &mf_real) {
    mf_real.copy(*mf_pointer_real,0,0,1);
    delete mf_pointer_real;
  }

  // If mf_imag is not defined on the particle_box_array, then we need                                                                                                                                   
  // to copy here from mf_pointer_imag into mf_imag. I believe that we don't                                                                                                                                  
  // need any information in ghost cells so we don't copy those.                                                                                                                                                 
  if (mf_pointer_imag != &mf_imag) {
    mf_imag.copy(*mf_pointer_imag,0,0,1);
    delete mf_pointer_imag;
  }

  if (m_verbose > 1) {
    Real stoptime = amrex::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "FDMParticleContainer::DepositFDMParticles time: " << stoptime << '\n';
  }
}

amrex::Real
FDMParticleContainer::estTimestepFDM(amrex::MultiFab&       phi,
				     amrex::Real              a,
				     int                    lev,
				     amrex::Real            cfl) const
{
  BL_PROFILE("FDMParticleContainer::estTimestep()");
  amrex::Real            dt               = 1e50;
  BL_ASSERT(lev >= 0);

  if (this->GetParticles().size() == 0)
    return dt;

  const amrex::Real      strttime         = amrex::ParallelDescriptor::second();
  const amrex::Geometry& geom             = this->m_gdb->Geom(lev);
  const ParticleLevel&   pmap             = this->GetParticles(lev);
  int                    tnum             = 1;

#ifdef _OPENMP
  tnum = omp_get_max_threads();
#endif

  amrex::Vector<amrex::Real> ldt(tnum,1e50);

  long num_particles_at_level = 0;
  amrex::MultiFab* phi_pointer;
  if (this->OnSameGrids(lev, phi))
    {
      phi_pointer = 0;
    }
  else
    {
      phi_pointer = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev),
				       this->m_gdb->ParticleDistributionMap(lev),
				       phi.nComp(), phi.nGrow());

      phi_pointer->copy(phi,0,0,1);
      phi_pointer->FillBoundary(geom.periodicity()); // DO WE NEED GHOST CELLS FILLED ???                                                                                                                         
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MyConstParIter pti(*this, lev); pti.isValid(); ++pti) {
    const int grid = pti.index();
    const AoS&       pbox = pti.GetArrayOfStructs();
    const int        n    = pbox.size();
    const amrex::FArrayBox& phifab = (phi_pointer) ? (*phi_pointer)[grid] : phi[grid];

    num_particles_at_level += n;
    for (int i = 0; i < n; i++) {
      const ParticleType& p = pbox[i];

      if (p.id() <= 0) continue;

      amrex::Real vel_square = (p.rdata(1)*p.rdata(1)+p.rdata(2)*p.rdata(2)+p.rdata(3)*p.rdata(3));
      amrex::Real dt_part = (vel_square > 0) ? (cfl * 2.0 * a * a / vel_square) : 1e50;
      amrex::IntVect cell = this->Index(p, lev);
      const amrex::Real pot = phifab(cell,0);
      if (pot > 0)
	dt_part = std::min( dt_part, cfl / pot );

      int tid = 0;

#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      ldt[tid] = std::min(dt_part, ldt[tid]);
    }
  }

  if (phi_pointer) delete phi_pointer;

  for (int i = 0; i < ldt.size(); i++)
    dt = std::min(dt, ldt[i]);

  amrex::ParallelDescriptor::ReduceRealMin(dt);
  //                                                                                                                                                                                                             
  // Set dt negative if there are no particles at this level.                                                                                                                                                    
  //                                                                                                                                                                                                             
  amrex::ParallelDescriptor::ReduceLongSum(num_particles_at_level);

  if (num_particles_at_level == 0) dt = -1.e50;
  if (this->m_verbose > 1)
    {
      amrex::Real stoptime = amrex::ParallelDescriptor::second() - strttime;

      amrex::ParallelDescriptor::ReduceRealMax(stoptime,amrex::ParallelDescriptor::IOProcessorNumber());

      if (amrex::ParallelDescriptor::IOProcessor())
        {
	  std::cout << "FDMParticleContainer::estTimestep() time: " << stoptime << '\n';
        }
    }

  return dt;
}


void
FDMParticleContainer::InitCosmo1ppcMultiLevel(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf,
					      amrex::Vector<amrex::MultiFab*>& phase,
					      const Real gamma_ax, const Real particleMass, 
					      BoxArray &baWhereNot, int lev, int nlevs)
{
    BL_PROFILE("FDMParticleContainer::InitCosmo1ppcMultiLevel()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(lev);
    const Real*     dx       = geom.CellSize();

    static Vector<int> calls;

    calls.resize(nlevs);

    calls[lev]++;

    if (calls[lev] > 1) return;

    Vector<ParticleLevel>& particles = this->GetParticles();

    particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    particles.resize(nlevs);

    ParticleType p;

    //
    // The mf should be initialized according to the ics...
    //
    int outside_counter=0;
    long outcount[3]={0,0,0};
    long outcountminus[3]={0,0,0};
    long totalcount=0;

    Real Amp;

    for (MFIter mfi(*(mf[lev])); mfi.isValid(); ++mfi)
    {
      FArrayBox&  myFab  = (*(mf[lev]))[mfi];
      FArrayBox&  phaseFab  = (*(phase[lev]))[mfi];
      const Box&  vbx    = mfi.validbox();
      const int  *fab_lo = vbx.loVect();
      const int  *fab_hi = vbx.hiVect();
      ParticleLocData pld;
      for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
	  for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
            	    IntVect indices(D_DECL(ix, jx, kx));
		    totalcount++;
		    if (baWhereNot.contains(indices)) 
		    {
                       continue;
		    }

	            for (int n = 0; n < BL_SPACEDIM; n++)
	            {
		      //
		      // Set positions (1 p per cell in the center of the cell)
		      //
		      p.pos(n) = geom.ProbLo(n) + 
			(indices[n]+Real(0.5))*dx[n];
		      //
		      // Set velocities
		      //
		      p.rdata(n+1) = myFab(indices,n+1);
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.rdata(0)  = particleMass;
	            p.id()      = ParticleType::NextID();
	            p.cpu()     = MyProc;

		    Amp = std::sqrt(myFab(indices,0));

		    //set phase
		    p.rdata( 4) = phaseFab(indices,0);
		    //set amplitude
		    p.rdata( 5) = pow(2.0*gamma_ax*Amp,1.5);
		    p.rdata( 6) = 0.0;
		    //set width
		    p.rdata( 7) = gamma_ax;
		    //set Jacobian qq
		    p.rdata( 8) = Amp;
		    p.rdata( 9) = 0.0;
		    p.rdata(10) = 0.0;
		    p.rdata(11) = 0.0;
		    p.rdata(12) = Amp;
		    p.rdata(13) = 0.0;
		    p.rdata(14) = 0.0;
		    p.rdata(15) = 0.0;
		    p.rdata(16) = Amp;
		    //set Jacobian pq
		    p.rdata(17) = 0.0;
		    p.rdata(18) = 0.0;
		    p.rdata(19) = 0.0;
		    p.rdata(20) = 0.0;
		    p.rdata(21) = 0.0;
		    p.rdata(22) = 0.0;
		    p.rdata(23) = 0.0;
		    p.rdata(24) = 0.0;
		    p.rdata(25) = 0.0;
		    //set Jacobian qp
		    p.rdata(26) = 0.0;
		    p.rdata(27) = 0.0;
		    p.rdata(28) = 0.0;
		    p.rdata(29) = 0.0;
		    p.rdata(30) = 0.0;
		    p.rdata(31) = 0.0;
		    p.rdata(32) = 0.0;
		    p.rdata(33) = 0.0;
		    p.rdata(34) = 0.0;
		    //set Jacobian pp
		    p.rdata(35) = Amp;
		    p.rdata(36) = 0.0;
		    p.rdata(37) = 0.0;
		    p.rdata(38) = 0.0;
		    p.rdata(39) = Amp;
		    p.rdata(40) = 0.0;
		    p.rdata(41) = 0.0;
		    p.rdata(42) = 0.0;
		    p.rdata(43) = Amp;


	
	            if (!this->Where(p, pld))
                    {
      		        this->PeriodicShift(p);

                        if (!this->Where(p, pld))
                            amrex::Abort("FDMParticleContainer::InitCosmo1ppcMultiLevel():invalid particle");
                    }

		    BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= m_gdb->finestLevel());
		    //handle particles that ran out of this level into a finer one. 
		    if (baWhereNot.contains(pld.m_cell))
		    {
		      outside_counter++;
		      ParticleType newp[8];
                      ParticleLocData new_pld;
		      for (int i=0;i<8;i++)
		      {
                          newp[i].rdata(0)   = particleMass/8.0;
                          newp[i].id()       = ParticleType::NextID();
                          newp[i].cpu()      = MyProc;
                          for (int dim=0;dim<BL_SPACEDIM;dim++)
                          {
                              newp[i].pos(dim)=p.pos(dim)+(2*((i/(1 << dim)) % 2)-1)*dx[dim]/4.0;
                              newp[i].rdata(dim+1)=p.rdata(dim+1);
                          }
                          
                          if (!this->Where(newp[i], new_pld))
                          {
                              this->PeriodicShift(newp[i]);
                              
                              if (!this->Where(newp[i], new_pld))
                                  amrex::Abort("FDMParticleContainer::InitCosmo1ppcMultiLevel():invalid particle");
                          }
                          particles[new_pld.m_lev][std::make_pair(new_pld.m_grid, 
                                                                  new_pld.m_tile)].push_back(newp[i]);
		      }
		      
		    }
	            
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
		    else
                        particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
                }
            }
        }
    }
    Redistribute();
}

void
FDMParticleContainer::InitCosmo1ppc(MultiFab& mf, const Real vel_fac[], const Real particleMass)
{
    BL_PROFILE("FDMParticleContainer::InitCosmo1ppc()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(0);
    const Real*     dx       = geom.CellSize();

    Vector<ParticleLevel>& particles = this->GetParticles();

    particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < particles.size(); lev++)
    {
        BL_ASSERT(particles[lev].empty());
    }

    ParticleType      p;
    ParticleLocData   pld;
    Real              disp[BL_SPACEDIM];
    const Real        len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                                  geom.ProbLength(1),
                                                  geom.ProbLength(2)) };
    //
    // The grid should be initialized according to the ics...
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
	const Box&  vbx    = mfi.validbox();
        const int  *fab_lo = vbx.loVect();
        const int  *fab_hi = vbx.hiVect();

        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
            for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
            	    IntVect indices(D_DECL(ix, jx, kx));

	            for (int n = 0; n < BL_SPACEDIM; n++)
	            {
                        disp[n] = myFab(indices,n);
                        //
			// Start with homogeneous distribution (for 1 p per cell in the center of the cell),
                        // then add the displacement (input values weighted by domain length).
                        //
	                p.pos(n) = geom.ProbLo(n) + 
                            (indices[n]+Real(0.5))*dx[n] +
                            disp[n] * len[n];
                        //
			// Set the velocities.
                        //
	                p.rdata(n+1) = disp[n] * vel_fac[n];
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.rdata(0)  = particleMass;
	            p.id()      = ParticleType::NextID();
	            p.cpu()     = MyProc;
	
	            if (!this->Where(p, pld))
                    {
      		        this->PeriodicShift(p);
                        
                        if (!this->Where(p, pld))
                            amrex::Abort("FDMParticleContainer::InitCosmo1ppc(): invalid particle");
		    }

	            BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= this->finestLevel());
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
	            particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
                }
            }
        }
    }
}

void
FDMParticleContainer::InitCosmo(
            MultiFab& mf, const Real vel_fac[], const Vector<int> n_part, const Real particleMass)
{
    Real shift[] = {0,0,0};
    InitCosmo(mf, vel_fac, n_part, particleMass, shift);
}

void
FDMParticleContainer::InitCosmo(
            MultiFab& mf, const Real vel_fac[], const Vector<int> n_part, const Real particleMass, const Real shift[])
{
    BL_PROFILE("FDMParticleContainer::InitCosmo()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(0);

    Vector<ParticleLevel>& particles = this->GetParticles();

    particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < particles.size(); lev++)
    {
        BL_ASSERT(particles[lev].empty());
    }

    const Real len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                           geom.ProbLength(1),
                                           geom.ProbLength(2)) };
    //
    // Print the grids as a sanity check.
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
	const Box&  vbx    = mfi.validbox();
        const int  *fab_lo = vbx.loVect();
        const int  *fab_hi = vbx.hiVect();
        if (vbx.isEmpty())
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           amrex::Error("Empty box in InitCosmo ");
        }
        if (!geom.Domain().contains(vbx))
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           amrex::Error("Box in InitCosmo not contained in domain");
        }
    }
    //
    // We will need one ghost cell, so check wether we have one.
    //
    if (mf.nGrow() < 1)
        amrex::Abort("FDMParticleContainer::InitCosmo: mf needs at least one correctly filled ghost zone!");

    if ( !(n_part[0] == n_part[1] && n_part[1] == n_part[2]) )
    {
	    std::cout << '\n' << '\n';
	    std::cout << "Your particle lattice will have different spacings in the spatial directions!" << '\n';
	    std::cout << "You might want to change the particle number or the algorithm... ;)" << '\n';
	    std::cout << '\n' << '\n';
    }
    //
    // Place the particles evenly spaced in the problem domain.
    // Not perfectly fast - but easy
    //
    Real         pos[BL_SPACEDIM];
    ParticleType p;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box&  box     = mfi.validbox();
        RealBox     gridloc = RealBox(box, geom.CellSize(), geom.ProbLo());
	const Real* xlo     = gridloc.lo();
	const Real* xhi     = gridloc.hi();

        ParticleLocData pld;
        for (int k = 0; k < n_part[2]; k++)
        {
	    for (int j = 0; j < n_part[1]; j++)
            {
	        for (int i = 0; i < n_part[0]; i++)
                {
		    bool    isInValidBox = true;
            	    IntVect indices(D_DECL(i, j, k));

		    for (int n = 0; n < BL_SPACEDIM; n++)
                    {
  		        pos[n] = geom.ProbLo(n)
		               + (indices[n] + Real(0.5))*len[n]/n_part[n]
			       + shift[n];
                        //
			// Make sure particle is not on a boundary...
                        //
			pos[n] += 1e-14 * (geom.ProbHi(n) - geom.ProbLo(n));

			isInValidBox = isInValidBox 
				     && (pos[n] > xlo[n]) 
				     && (pos[n] < xhi[n]);
		    }

		    if (isInValidBox)
                    {
                        D_TERM(p.pos(0) = pos[0];,
                               p.pos(1) = pos[1];,
                               p.pos(2) = pos[2];);
                        //
		        // Set the mass of the particle.
                        //
	                p.rdata(0)  = particleMass;
	                p.id()      = ParticleType::NextID();
	                p.cpu()     = MyProc;

	                if (!this->Where(p, pld))
                        {
      		            this->PeriodicShift(p);

                            if (!this->Where(p, pld))
                                amrex::Abort("FDMParticleContainer::InitCosmo(): invalid particle");
		        }

	                BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= m_gdb->finestLevel());
	                //
	                // Add it to the appropriate PBox at the appropriate level.
	                //
	                particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
		    }
	        }
	    }
        }
    }

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Done with equidistant placement" << '\n';
    }
    //
    // Let Redistribute() sort out where the particles belong.
    //
    Redistribute();

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Redistribute done" << '\n';
    }

    BL_ASSERT(OK());
    //
    // FIXME: Will we ever need initial particles in grids deeper than 0?!
    //
    ParticleLevel& pmap = this->GetParticles(0);
    //
    // Make sure, that mf and m_gdb.boxArray(0) are defined on the same boxarray.
    //
    ParticleLocData pld;
     for (auto& kv : pmap) {
        const int        grid    = kv.first.first;
        AoS&             pbox    = kv.second.GetArrayOfStructs();
        const int        n       = pbox.size();
        const FArrayBox& dfab    = mf[grid];

#ifdef _OPENMP
#pragma omp parallel for private(pld)
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.id() <= 0) continue;

            Real disp[BL_SPACEDIM];
            //
	    // Do CIC interpolation onto the particle positions.
	    // For CIC we need one ghost cell!
            //
            ParticleType::GetGravity(dfab, m_gdb->Geom(0), p, disp);

            D_TERM(p.pos(0) += len[0]*disp[0];,
                   p.pos(1) += len[1]*disp[1];,
                   p.pos(2) += len[2]*disp[2];);
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            D_TERM(p.rdata(1) = vel_fac[0]*disp[0];,
                   p.rdata(2) = vel_fac[1]*disp[1];,
                   p.rdata(3) = vel_fac[2]*disp[2];);

            if (!this->Where(p, pld))
            {
	        this->PeriodicShift(p);

                if (!this->Where(p, pld))
                    amrex::Abort("FDMParticleContainer::InitCosmo(): invalid particle");
	    }

            this->Reset(p, true);
        }
    }
    //
    // Let Redistribute() sort out where the particles now belong.
    //
    Redistribute();

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Done with particle displacement" << '\n';
    }

    if (m_verbose > 1)
    {
        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "InitCosmo() done time: " << runtime << '\n';
        }
    }
}

void 
FDMParticleContainer::InitFromBinaryMortonFile(const std::string& particle_directory,
						      int nextra, int skip_factor) {
  BL_PROFILE("FDMParticleContainer::InitFromBinaryMortonFile");
  
  ParticleMortonFileHeader hdr;
  ReadHeader(particle_directory, "Header", hdr);    
  
  uint64_t num_parts = hdr.NP;
  int DM             = hdr.DM;
  int NX             = hdr.NX;
  int float_size     = hdr.SZ;
  int num_files      = hdr.NF;
  size_t psize       = (DM + NX) * float_size;
  
  std::string particle_file_base = particle_directory + "/particles.";
  std::vector<std::string> file_names;
  for (int i = 0; i < num_files; ++i)
    file_names.push_back(get_file_name(particle_file_base, i));
  
  const int lev = 0;
  const BoxArray& ba = ParticleBoxArray(lev);
  int num_boxes = ba.size();
  uint64_t num_parts_per_box  = num_parts / num_boxes;
  uint64_t num_parts_per_file = num_parts / num_files;
  uint64_t num_bytes_per_file = num_parts_per_file * psize;
  
  std::vector<BoxMortonKey> box_morton_keys(num_boxes);
  for (int i = 0; i < num_boxes; ++i) {
    const Box& box = ba[i];
    unsigned int x = box.smallEnd(0);
    unsigned int y = box.smallEnd(1);
    unsigned int z = box.smallEnd(2);
    box_morton_keys[i].morton_id = get_morton_index(x, y, z);
    box_morton_keys[i].box_id = i;
  }
  
  std::sort(box_morton_keys.begin(), box_morton_keys.end(), by_morton_id());
  
  std::vector<int> file_indices(num_boxes);
  for (int i = 0; i < num_boxes; ++i)
    file_indices[box_morton_keys[i].box_id] = i;
  
  ParticleType p;
  for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {  // no tiling
    Box tile_box = mfi.tilebox();      
    const int grid = mfi.index();
    const int tile = mfi.LocalTileIndex();      
    auto& particles = GetParticles(lev);
    
    uint64_t start    = file_indices[grid]*num_parts_per_box;
    uint64_t stop     = start + num_parts_per_box;

    int file_num      = start / num_parts_per_file;
    uint64_t seek_pos = (start * psize ) % num_bytes_per_file;
    std::string file_name = file_names[file_num];
    
    std::ifstream ifs;
    ifs.open(file_name.c_str(), std::ios::in|std::ios::binary);
    if ( not ifs ) {
      amrex::Print() << "Failed to open file " << file_name << " for reading. \n";
      amrex::Abort();
    } 

    ifs.seekg(seek_pos, std::ios::beg);
    
    for (uint64_t i = start; i < stop; ++i) {
      int next_file = i / num_parts_per_file;
      if (next_file != file_num) {
	file_num = next_file;
	file_name = file_names[file_num];
	ifs.close();
	ifs.open(file_name.c_str(), std::ios::in|std::ios::binary);
	if ( not ifs ) {
	  amrex::Print() << "Failed to open file " << file_name << " for reading. \n";
	  amrex::Abort();
	}
      }

      float fpos[DM];
      float fextra[NX];
      ifs.read((char*)&fpos[0],   DM*sizeof(float));
      ifs.read((char*)&fextra[0], NX*sizeof(float));
      
      if ( (i - start) % skip_factor == 0 ) {
	AMREX_D_TERM(p.m_rdata.pos[0] = fpos[0];,
		     p.m_rdata.pos[1] = fpos[1];,
		     p.m_rdata.pos[2] = fpos[2];);
	
	for (int comp = 0; comp < NX; comp++)
	  p.m_rdata.arr[BL_SPACEDIM+comp] = fextra[comp];
	
	p.m_rdata.arr[BL_SPACEDIM] *= skip_factor;
	
	p.m_idata.id  = ParticleType::NextID();
	p.m_idata.cpu = ParallelDescriptor::MyProc();
	particles[std::make_pair(grid, tile)].push_back(p);
      }
    }    
  }
  
  Redistribute();
}

void
FDMParticleContainer::InitVarCount (MultiFab& mf, long num_particle_fdm, BoxArray &baWhereNot, int lev, int nlevs)
{
  const int       MyProc      = ParallelDescriptor::MyProc();
  const Geometry& geom        = m_gdb->Geom(lev);
  const Real*     dx          = geom.CellSize();

  static Vector<int> calls;
  calls.resize(nlevs);
  calls[lev]++;
  if (calls[lev] > 1) return;
  Vector<ParticleLevel>& particles = this->GetParticles();
  int npart;
  Real r;
  Real factor = mf.norm1(0,0)/num_particle_fdm; //compute density per particle                                                                                             
  particles.reserve(15);  // So we don't ever have to do any copying on a resize.                                                                                                                        
  particles.resize(nlevs);

  for (int i = 0; i < particles.size(); i++)
    {
      BL_ASSERT(particles[i].empty());
    }

  ParticleType p;

  amrex::InitRandom(MyProc);
  //                                                                                                                                                                                                             
  // The grid should be initialized according to the ics...                                                                                                                                                      
  //                                                                                                                                                                                                             
  for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {

      FArrayBox&  myFab  = mf[mfi];
      const Box&  vbx    = mfi.validbox();
      const int  *fab_lo = vbx.loVect();
      const int  *fab_hi = vbx.hiVect();
      ParticleLocData pld;

      for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
	  for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
	      for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
		  IntVect indices(D_DECL(ix, jx, kx));

		  //compute the floor of number of particles in the cell                                                                                                                                         
		  npart = myFab(indices,0)/factor;

		  //roll dice to decide whether to add another particle -----                                                                                                                                    
		  r = amrex::Random();

		  if(r < myFab(indices,0)/factor-npart){
		    npart++;
		  }
		  // --------------------------------------------------------                                                                                                                                    

		  for (int ipart = 0; ipart < npart; ipart++){



		    for (int n = 0; n < BL_SPACEDIM; n++)
		      {
                                                do
						  {
						    r = amrex::Random();
						  }
                                                while (r == 0 || r == 1);

                                                p.pos(n) = geom.ProbLo(n) +
						  (indices[n]+r)*dx[n];
		      }

		    // set mass                                                                                                                                                                
		    p.rdata(0) = factor * dx[0] * dx[1] * dx[2];
		    // set velocity                                                                                                                                                            
		    for (int n = 0; n < BL_SPACEDIM; n++)
		      {
			p.rdata(n+1) = myFab(indices, n+1);
		      }

		    p.id()      = ParticleType::NextID();
		    p.cpu()     = MyProc;

		    if (!this->Where(p,pld))
		      {
			this->PeriodicShift(p);

			if (!this->Where(p,pld))
			  amrex::Abort("ParticleContainer<N>::InitVarCount(): invalid particle");
		      }

		    BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= m_gdb->finestLevel());
                    //handle particles that ran out of this level into a finer one.                                                                                                                                
                    if (baWhereNot.contains(pld.m_cell))
		      {
			ParticleLocData new_pld;
			if (!this->Where(p, new_pld))
			  {
			    this->PeriodicShift(p);
			    
			    if (!this->Where(p, new_pld))
			      amrex::Abort("FDMParticleContainer::InitVarCount():invalid particle");
			  }
			particles[new_pld.m_lev][std::make_pair(new_pld.m_grid,
								new_pld.m_tile)].push_back(p);
			
		      }


		    //                                                                                                                                                                         
		    // Add it to the appropriate PBox at the appropriate level.                                                                                                                
		    //                                                                                                                                                                         
		    else
		      particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
		  }
                }
            }
        }
    }

  if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
      std::cout << "Done with Gaussian Beam initilization" << '\n';
    }
  //
  // Let Redistribute() sort out where the particles belong.
  //
  Redistribute();
  
  if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
      std::cout << "Redistribute done" << '\n';
    }
}

void
FDMParticleContainer::InitGaussianBeams (long num_particle_fdm, int lev, int nlevs, const Real hbaroverm, const Real sigma_ax, const Real gamma_ax, const Real fact, const Real alpha, const Real a)
{
  const int       MyProc      = ParallelDescriptor::MyProc();
  const int       nprocs      = ParallelDescriptor::NProcs();
  const Geometry& geom        = m_gdb->Geom(lev);
  const Real*     dx          = geom.CellSize();

  static Vector<int> calls;
  calls.resize(nlevs);
  calls[lev]++;
  if (calls[lev] > 1) return;
  Vector<ParticleLevel>& particles = this->GetParticles();

  int  npart     = num_particle_fdm;
  int  npart_tot = nprocs*npart; //Each processor initializes num_particle_fdm beams
  Real q[]  = {(geom.ProbHi(0)+geom.ProbLo(0))/2.0, (geom.ProbHi(1)+geom.ProbLo(1))/2.0, (geom.ProbHi(2)+geom.ProbLo(2))/2.0};
  Real p[]  = {0.0,0.0,0.0};
  Real q0[]  = {(geom.ProbHi(0)+geom.ProbLo(0))/2.0, (geom.ProbHi(1)+geom.ProbLo(1))/2.0, (geom.ProbHi(2)+geom.ProbLo(2))/2.0};
  Real p0[] = {0.0,0.0,0.0};
  Real r, theta, phi, Amp;

  particles.reserve(15);  // So we don't ever have to do any copying on a resize.
  particles.resize(nlevs);

  for (int i = 0; i < particles.size(); i++)
    {
      BL_ASSERT(particles[i].empty());
    }

  ParticleType part;
  ParticleLocData pld;

  amrex::InitRandom(MyProc);

  for(int index=0;index<npart;index++){
	
    r = sqrt(-log(1-amrex::Random()));
    theta  = 2.0*M_PI*amrex::Random();
    q[0] = r*cos(theta)*sqrt((alpha+gamma_ax)/gamma_ax/alpha) + q0[0];
    p[0] = r*sin(theta)*2.0*sqrt(alpha+gamma_ax)*hbaroverm + p0[0];
	
    r = sqrt(-log(1-amrex::Random()));
    theta  = 2.0*M_PI*amrex::Random();
    q[1] = r*cos(theta)*sqrt((alpha+gamma_ax)/gamma_ax/alpha) + q0[1];
    p[1] = r*sin(theta)*2.0*sqrt(alpha+gamma_ax)*hbaroverm + p0[1];
	
    r = sqrt(-log(1-amrex::Random()));
    theta  = 2.0*M_PI*amrex::Random();
    q[2] = r*cos(theta)*sqrt((alpha+gamma_ax)/gamma_ax/alpha) + q0[2];
    p[2] = r*sin(theta)*2.0*sqrt(alpha+gamma_ax)*hbaroverm + p0[2];
	
    phi  = ( (p[0]*alpha+p0[0]*gamma_ax)*(q[0]-q0[0]) + (p[1]*alpha+p0[1]*gamma_ax)*(q[1]-q0[1]) + (p[2]*alpha+p0[2]*gamma_ax)*(q[2]-q0[2]) )/(alpha+gamma_ax);
    Amp  = 2.0*(alpha+gamma_ax)/sqrt(alpha*gamma_ax)/M_PI/sqrt(2*gamma_ax/M_PI)/pow(npart_tot,2.0/3.0);
    Amp /= sqrt(2.0*alpha/M_PI);
    Amp *= pow(100.0*fact,1.0/3.0);

    if(q[0]>geom.ProbLo(0) && q[0]<geom.ProbHi(0) && q[1]>geom.ProbLo(1) && q[1]<geom.ProbHi(1) && q[2]>geom.ProbLo(2) && q[2]<geom.ProbHi(2)){

      part.id()      = ParticleType::NextID();
      part.cpu()     = MyProc;

      //set position
      for (int n = 0; n < BL_SPACEDIM; n++)
	part.pos( n) = q[n];
      //set mass                                                                                                                                                                
      part.rdata( 0) =  1.0/(npart_tot * dx[0] * dx[1] * dx[2]);
      //set velocity
      part.rdata( 1) = p[0]/a;
      part.rdata( 2) = p[1]/a;
      part.rdata( 3) = p[2]/a;
      //set phase
      part.rdata( 4) = phi;
      //set amplitude
      part.rdata( 5) = pow(2.0*gamma_ax*Amp,1.5);
      part.rdata( 6) = 0.0;
      //set width
      part.rdata( 7) = gamma_ax;
      //set Jacobian qq
      part.rdata( 8) = Amp;
      part.rdata( 9) = 0.0;
      part.rdata(10) = 0.0;
      part.rdata(11) = 0.0;
      part.rdata(12) = Amp;
      part.rdata(13) = 0.0;
      part.rdata(14) = 0.0;
      part.rdata(15) = 0.0;
      part.rdata(16) = Amp;
      //set Jacobian pq
      part.rdata(17) = 0.0;
      part.rdata(18) = 0.0;
      part.rdata(19) = 0.0;
      part.rdata(20) = 0.0;
      part.rdata(21) = 0.0;
      part.rdata(22) = 0.0;
      part.rdata(23) = 0.0;
      part.rdata(24) = 0.0;
      part.rdata(25) = 0.0;
      //set Jacobian qp
      part.rdata(26) = 0.0;
      part.rdata(27) = 0.0;
      part.rdata(28) = 0.0;
      part.rdata(29) = 0.0;
      part.rdata(30) = 0.0;
      part.rdata(31) = 0.0;
      part.rdata(32) = 0.0;
      part.rdata(33) = 0.0;
      part.rdata(34) = 0.0;
      //set Jacobian pp
      part.rdata(35) = Amp;
      part.rdata(36) = 0.0;
      part.rdata(37) = 0.0;
      part.rdata(38) = 0.0;
      part.rdata(39) = Amp;
      part.rdata(40) = 0.0;
      part.rdata(41) = 0.0;
      part.rdata(42) = 0.0;
      part.rdata(43) = Amp;

      if (!this->Where(part,pld))
	amrex::Abort("ParticleContainer<N>::InitGaussianBeams(): invalid particle");

      //add particle
      particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(part);
    }
    else
      index--;
  }

  if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
      std::cout << "Done with Gaussian Beam initilization" << '\n';
    }
  //
  // Let Redistribute() sort out where the particles belong.
  //
  Redistribute();
  
  if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
      std::cout << "Redistribute done" << '\n';
    }
  
}

amrex::Real
FDMParticleContainer::generateGaussianNoise(const amrex::Real &mean, const amrex::Real &stdDev) {

  static bool hasSpare = false;
  static amrex::Real spare;

  if(hasSpare) {
    hasSpare = false;
    return mean + stdDev * spare;
  }

  hasSpare = true;
  static amrex::Real u, v, s;
  do {
    u = (rand() / ((amrex::Real) RAND_MAX)) * 2.0 - 1.0;
    v = (rand() / ((amrex::Real) RAND_MAX)) * 2.0 - 1.0;
    s = u * u + v * v;
  }
  while( (s >= 1.0) || (s == 0.0) );

  s = sqrt(-2.0 * log(s) / s);
  spare = v * s;
  return mean + stdDev * u * s;
}

#endif /*FDM*/
