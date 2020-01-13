#include <stdint.h>

#include "DarkMatterParticleContainer.H"
#include "dm_F.H"

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
DarkMatterParticleContainer::moveKickDrift (amrex::MultiFab&       acceleration,
		                            int                    lev,
                    			    amrex::Real            dt,
		                	    amrex::Real            a_old,
					    amrex::Real            a_half,
					    int                    where_width)
{
    BL_PROFILE("DarkMatterParticleContainer::moveKickDrift()");

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

           update_dm_particles(&Np, particles.data(),
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
                        int grid = kv.first.first;
                        
                        
                        std::cout << "Oops -- removing particle " << p << " " << this->Index(p, lev) << " " << lev << " " << (this->m_gdb->ParticleBoxArray(lev))[grid] << " " << where_width << std::endl;
                        amrex::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
        }
    }
}

void
DarkMatterParticleContainer::moveKick (MultiFab&       acceleration,
                                       int             lev,
                                       Real            dt,
                                       Real            a_new,
                                       Real            a_half) 
{
    BL_PROFILE("DarkMatterParticleContainer::moveKick()");

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

           update_dm_particles(&Np, particles.data(),
                               (*ac_ptr)[pti].dataPtr(),
                               ac_box.loVect(), ac_box.hiVect(),
                               plo,dx,dt,a_half,a_new,&do_move);
        }
    }
    
    if (ac_ptr != &acceleration) delete ac_ptr;
}

void
DarkMatterParticleContainer::InitCosmo1ppcMultiLevel(
                        MultiFab& mf, const Real disp_fac[], const Real vel_fac[], 
                        const Real particleMass, int disp_idx, int vel_idx, 
                        BoxArray &baWhereNot, int lev, int nlevs)
{
    BL_PROFILE("DarkMatterParticleContainer::InitCosmo1ppcMultiLevel()");
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
    Real         disp[BL_SPACEDIM];
    Real         vel[BL_SPACEDIM];
    
    Real 	mean_disp[BL_SPACEDIM]={D_DECL(0,0,0)};


    //
    // The mf should be initialized according to the ics...
    //
    int outside_counter=0;
    long outcount[3]={0,0,0};
    long outcountminus[3]={0,0,0};
    long totalcount=0;
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
    		    totalcount++;
    		    if (baWhereNot.contains(indices)) 
    		    {
                       continue;
    		    }

    	            for (int n = 0; n < BL_SPACEDIM; n++)
    	            {
                        disp[n] = myFab(indices,disp_idx+n);
                        //
    			// Start with homogeneous distribution (for 1 p per cell in the center of the cell),
    			//
    	                p.pos(n) = geom.ProbLo(n) + 
                            (indices[n]+Real(0.5))*dx[n];
    			if(disp[n]*disp_fac[n]>dx[n]/2.0)
    			  outcount[n]++;
    			if(disp[n]*disp_fac[n]<-dx[n]/2.0)
    			  outcountminus[n]++;
    			mean_disp[n]+=fabs(disp[n]);
    			//
                        // then add the displacement (input values weighted by domain length).
                        //
    	                p.pos(n) += disp[n] * disp_fac[n];

                        //
    			// Set the velocities.
                        //
                        vel[n] = myFab(indices,vel_idx+n);
    	                p.rdata(n+1) = vel[n] * vel_fac[n];
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
                            amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppcMultiLevel():invalid particle");
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
                                  amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppcMultiLevel():invalid particle");
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
DarkMatterParticleContainer::InitCosmo1ppc(MultiFab& mf, const Real vel_fac[], const Real particleMass)
{
    BL_PROFILE("DarkMatterParticleContainer::InitCosmo1ppc()");
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
                            amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppc(): invalid particle");
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
DarkMatterParticleContainer::InitCosmo(
            MultiFab& mf, const Real vel_fac[], const Vector<int> n_part, const Real particleMass)
{
    Real shift[] = {0,0,0};
    InitCosmo(mf, vel_fac, n_part, particleMass, shift);
}

void
DarkMatterParticleContainer::InitCosmo(
            MultiFab& mf, const Real vel_fac[], const Vector<int> n_part, const Real particleMass, const Real shift[])
{
    BL_PROFILE("DarkMatterParticleContainer::InitCosmo()");
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
        amrex::Abort("DarkMatterParticleContainer::InitCosmo: mf needs at least one correctly filled ghost zone!");

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
                                amrex::Abort("DarkMatterParticleContainer::InitCosmo(): invalid particle");
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
                    amrex::Abort("DarkMatterParticleContainer::InitCosmo(): invalid particle");
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

/*
  Particle deposition
*/

void
DarkMatterParticleContainer::AssignDensityAndVels (Vector<std::unique_ptr<MultiFab> >& mf, int lev_min) const
{
     AssignDensity(mf, lev_min, BL_SPACEDIM+1);
}

void 
DarkMatterParticleContainer::InitFromBinaryMortonFile(const std::string& particle_directory,
						      int nextra, int skip_factor) {
  BL_PROFILE("DarkMatterParticleContainer::InitFromBinaryMortonFile");
  
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
#ifdef FDM
void
DarkMatterParticleContainer::InitGaussianBeams (long num_particle_dm, int lev, int nlevs, const Real fact, const Real alpha, const Real a)
{

  const int       MyProc      = ParallelDescriptor::MyProc();
  const int       nprocs      = ParallelDescriptor::NProcs();
  const Geometry& geom        = m_gdb->Geom(lev);

  static Vector<int> calls;
  calls.resize(nlevs);
  calls[lev]++;
  if (calls[lev] > 1) return;
  Vector<ParticleLevel>& particles = this->GetParticles();

  int  npart = num_particle_dm;
  int  npart_tot = nprocs*npart; //Each processor initializes num_particle_fdm beams
  Real q0[]  = {(geom.ProbHi(0)+geom.ProbLo(0))/2.0, (geom.ProbHi(1)+geom.ProbLo(1))/2.0, (geom.ProbHi(2)+geom.ProbLo(2))/2.0};
  Real q[]  = {(geom.ProbHi(0)+geom.ProbLo(0))/2.0, (geom.ProbHi(1)+geom.ProbLo(1))/2.0, (geom.ProbHi(2)+geom.ProbLo(2))/2.0};
  Real p0[] = {0.0,0.0,0.0};
  // Real sigma = 0.5/sqrt(alpha);/*remember that we square amplitude alpha->2*alpha*/
  Real sigma = (geom.ProbHi(0)+geom.ProbLo(0))/16.0;

  //calculate dm particle mass
  Real mass = 1.0/npart_tot;
  // mass *= pow(2.0*alpha/M_PI,-1.5);
  mass *= pow(1.0/sigma/sigma/2.0/M_PI,-1.5);
  mass *= 0.1*fact;

  // std::cout<<"meandens cdm = "<<fact<<std::endl;

  Real mass2 = fact*pow(geom.CellSize(1),3);

  mass*=0.5;
  mass2*=0.5;

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

    q[0] = generateGaussianNoise(q0[0],sigma);
    q[1] = generateGaussianNoise(q0[1],sigma);
    q[2] = generateGaussianNoise(q0[2],sigma);

    if(q[0]>geom.ProbLo(0) && q[0]<geom.ProbHi(0) && q[1]>geom.ProbLo(1) && q[1]<geom.ProbHi(1) && q[2]>geom.ProbLo(2) && q[2]<geom.ProbHi(2)){
      
      part.id()      = ParticleType::NextID();
      part.cpu()     = MyProc;
      
      //set position
      for (int n = 0; n < BL_SPACEDIM; n++)
  	part.pos( n) = q[n];
      //set mass
      part.rdata( 0) =  mass;
      //set velocity
      part.rdata( 1) = p0[0]/a;
      part.rdata( 2) = p0[1]/a;
      part.rdata( 3) = p0[2]/a;
      
      if (!this->Where(part,pld))
  	amrex::Abort("ParticleContainer<N>::InitGaussianBeams(): invalid particle");
      
      //add particle                                                                                                                                                                                            
      particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(part);
    }    
    else
      index--;
  }


    // for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    // {
    //     FArrayBox&  myFab  = mf[mfi];
    // 	const Box&  vbx    = mfi.validbox();
    //     const int  *fab_lo = vbx.loVect();
    //     const int  *fab_hi = vbx.hiVect();

    //     for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
    //     {
    //         for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
    //         {
    //             for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
    //             {
    //         	    IntVect indices(D_DECL(ix, jx, kx));

    // 	            for (int n = 0; n < BL_SPACEDIM; n++)
    // 	            {
    //                     // disp[n] = myFab(indices,n);
    //                     //
    // 			// Start with homogeneous distribution (for 1 p per cell in the center of the cell),
    //                     // then add the displacement (input values weighted by domain length).
    //                     //
    // 	                part.pos(n) = geom.ProbLo(n) + 
    // 			  (indices[n]+Real(0.5))*dx[n];// +
    // 			//                            disp[n] * len[n];
    //                     //
    // 			// Set the velocities.
    //                     //
    // 	                part.rdata(n+1) = 0.0;//disp[n] * vel_fac[n];
    // 	            }
    //                 //
    // 		    // Set the mass of the particle from the input value.
    //                 //
    // 	            part.rdata(0)  = mass2;
    // 	            part.id()      = ParticleType::NextID();
    // 	            part.cpu()     = MyProc;
	
    // 	            if (!this->Where(part, pld))
    //                 {
    //   		        this->PeriodicShift(part);
                        
    //                     if (!this->Where(part, pld))
    //                         amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppc(): invalid particle");
    // 		    }

    // 	            BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= this->finestLevel());
    // 	            //
    // 	            // Add it to the appropriate PBox at the appropriate level.
    // 	            //
    // 	            particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(part);
    //             }
    //         }
    //     }
    // }




  const Box& box = geom.Domain();
  for (int kx = box.smallEnd(2); kx <= box.bigEnd(2); kx++)
    for (int jx = box.smallEnd(1); jx <= box.bigEnd(1); jx++)
      for (int ix = box.smallEnd(0); ix <= box.bigEnd(0); ix++)
  	{
  	  IntVect indices(D_DECL(ix, jx, kx));
  	  // Real r[] = {(geom.ProbHi(0)-geom.ProbLo(0))/2.0-(geom.ProbLo(0) + (indices[0]+Real(0.5))*geom.CellSize(0)),
  	  // 	      (geom.ProbHi(1)-geom.ProbLo(1))/2.0-(geom.ProbLo(1) + (indices[1]+Real(0.5))*geom.CellSize(1)),
  	  // 	      (geom.ProbHi(2)-geom.ProbLo(2))/2.0-(geom.ProbLo(2) + (indices[2]+Real(0.5))*geom.CellSize(2))};
  	  // Real rsq = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];

  	  for (int n = 0; n < BL_SPACEDIM; n++)
  	    {
	      
  	      part.pos(n) = geom.ProbLo(n) + (indices[n]+Real(0.5))*geom.CellSize(n);
  	      // part.pos(n) += 0.5*r[n]*rsq*std::exp(-rsq/0.1); 
  	      part.rdata(n+1) = 0.0;
  	    }
  	  part.rdata(0)  = mass2;
  	  part.id()      = ParticleType::NextID();
  	  part.cpu()     = MyProc;

  	  if (!this->Where(part,pld))
  	    amrex::Abort("ParticleContainer<N>::InitGaussianBeams(): invalid particle");
	  
  	  //add particle                                                                                                                                                                                       
  	  particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(part);

  	}

  if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
      std::cout << "Done DM particle initilization for Gaussian Beam potential.\n";
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
DarkMatterParticleContainer::generateGaussianNoise(const amrex::Real &mean, const amrex::Real &stdDev) {

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


void
DarkMatterParticleContainer::InitSphericalCollapse (amrex::MultiFab& mf, int lev, int nlevs, int comp, Real ratio)
{

  /*Pure FDM, so no nbody particles needed*/
  if(ratio==1.0) return;

  const int       MyProc      = ParallelDescriptor::MyProc();
  const int       nprocs      = ParallelDescriptor::NProcs();
  const Geometry& geom        = m_gdb->Geom(lev);

  static Vector<int> calls;
  calls.resize(nlevs);
  calls[lev]++;
  if (calls[lev] > 1) return;
  Vector<ParticleLevel>& particles = this->GetParticles();

  const Real* dx = geom.CellSize();

  particles.reserve(15);  // So we don't ever have to do any copying on a resize.                                                                                                                           
  particles.resize(nlevs);

  for (int i = 0; i < particles.size(); i++)
    {
      BL_ASSERT(particles[i].empty());
    }

  const Real rc = pow(0.125*(geom.ProbHi(0)-geom.ProbLo(0)),2)/2.0;
  const Real center = 0.5*(geom.ProbHi(0)+geom.ProbLo(0));

  // Real comoving_OmM,comoving_OmB,comoving_h,Gconst;
  // fort_get_omm(&comoving_OmM );                                                                                                                                                                               
  // fort_get_omb(&comoving_OmB );                                                                                                                                                                               
  // fort_get_hubble(&comoving_h);                                                                                                                                                                               
  // fort_get_grav_const(&Gconst); 
  // const Real meandens = 3*comoving_h*100*comoving_h*100*comoving_OmM / (8*M_PI*Gconst);

  const Real meandens = 2.775e+11*0.674*0.674*0.315;
  // Real r[] = {0.0,0.0,0.0};
  // Real rsq,disp;
  Real r;
  ParticleType part;
  ParticleLocData pld;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
    	const Box&  vbx    = mfi.validbox();
        const int  *fab_lo = vbx.loVect();
        const int  *fab_hi = vbx.hiVect();

        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
	  for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
	    for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
	      {
		IntVect indices(D_DECL(ix, jx, kx));

		// rsq = 0.0;
		// for (int n = 0; n < BL_SPACEDIM; n++)
		//   {
		//     r[n]= (indices[n]+0.5)*dx[n]-0.5*(geom.ProbHi(n)+geom.ProbLo(n));
		//     rsq += r[n]*r[n];
		//   }
		// // disp = 0.1*(0.25*sqrt(M_PI*rc/rsq)*erf(sqrt(rsq/rc))-0.5*exp(-rsq/rc));
		// disp = 0.1*0.5*sqrt(M_PI*rc/rsq)*erf(sqrt(rsq/rc));

		for (int n = 0; n < BL_SPACEDIM; n++)
		  {
		    part.pos(n) = geom.ProbLo(n) + (indices[n]+Real(0.5))*dx[n];//-disp*r[n];
 		    // part.pos(n)+= (rand() / ((amrex::Real) RAND_MAX) - 0.5)*dx[n]/10000.0;
		    part.rdata(n+1) = 0.0;
		  }
                    //
    		    // Set the mass of the particle from the input value.
                    //
		if(ratio>=0.0 && ratio <=1.0){

		  // part.rdata(0)  = 0.0;
		  // for(int i=0;i<10;i++)
		  //   for(int j=0;j<10;j++)
		  //     for(int k=0;k<10;k++){
			// r = pow((indices[0]+Real(i)/100.0)*dx[0]-center,2)+pow((indices[1]+Real(j)/100.0)*dx[1]-center,2)+pow((indices[2]+Real(k)/100.0)*dx[2]-center,2);
		  r = pow((indices[0]+0.5)*dx[0]-center,2)+pow((indices[1]+0.5)*dx[1]-center,2)+pow((indices[2]+0.5)*dx[2]-center,2);
		  // r = pow((indices[0]+0.0+(rand() / ((amrex::Real) RAND_MAX)/1.0))*dx[0]-center,2)
		  //    +pow((indices[1]+0.0+(rand() / ((amrex::Real) RAND_MAX)/1.0))*dx[1]-center,2)
		  //    +pow((indices[2]+0.0+(rand() / ((amrex::Real) RAND_MAX)/1.0))*dx[2]-center,2);
		  // r = pow((indices[0]+0.0)*dx[0]-center,2)+pow((indices[1]+0.0)*dx[1]-center,2)+pow((indices[2]+0.0)*dx[2]-center,2);
			// part.rdata(0)  += (1.0-ratio)*(0.001*meandens*exp(-r/rc)+meandens)*dx[0]*dx[1]*dx[2]/pow(10.0,3);
		  part.rdata(0)  = (1.0-ratio)*(0.1*meandens*exp(-r/rc)+meandens)*dx[0]*dx[1]*dx[2];
		  // }
		  // part.rdata(0)  = (1-ratio)/ratio*myFab(indices,comp)*dx[0]*dx[1]*dx[2];
		  // part.rdata(0)  = (1.0-ratio)*meandens*dx[0]*dx[1]*dx[2];
		}else
		  amrex::Abort("ParticleContainer<N>::InitSphericalCollapse(): ratio needs to be between 0 and 1.");
		part.id()      = ParticleType::NextID();
		part.cpu()     = MyProc;
		
		if (!this->Where(part, pld))
		  {
		    this->PeriodicShift(part);
                    
		    if (!this->Where(part, pld))
		      amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppc(): invalid particle");
		  }
		
		BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= this->finestLevel());
		//
		// Add it to the appropriate PBox at the appropriate level.
		//
		particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(part);
	      }
    }
    
    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
      std::cout << "Done DM particle initilization for Gaussian Beam potential.\n";
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


#endif
