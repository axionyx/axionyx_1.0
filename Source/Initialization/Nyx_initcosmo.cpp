#include <iomanip>

#ifdef GRAVITY
#include <Gravity.H>
extern "C"
{ void fort_get_grav_const(amrex::Real* Gconst); }
#endif

#include <Nyx.H>
#include <Nyx_F.H>
#include <NyxParticleContainer.H>
#ifdef FDM
#include "fdm.H"
#endif

using namespace amrex;

#ifdef GRAVITY
void Nyx::icReadAndPrepareFab(std::string mfDirName, int nghost, MultiFab &mf)
{
    if (level > 0 && nghost > 0)
    {
       std::cout << "Are sure you want to do what you are doing?" << std::endl;
       amrex::Abort();
    }

    //
    // Read from the directory into mf
    //
    MultiFab mf_read;

    if (!mfDirName.empty() && mfDirName[mfDirName.length()-1] != '/')
       mfDirName += '/';
    std::string Level = amrex::Concatenate("Level_", level, 1);
    mfDirName.append(Level);
    mfDirName.append("/Cell");

    VisMF::Read(mf_read,mfDirName.c_str());

    if (ParallelDescriptor::IOProcessor())
        std::cout << "mf read" << '\n';

    if (mf_read.contains_nan())
    {
        for (int i = 0; i < mf_read.nComp(); i++)
        {
            if (mf_read.contains_nan(i, 1))
            {
                std::cout << "Found NaNs in read_mf in component " << i << ". " << std::endl;
                amrex::Abort("Nyx::init_particles: Your initial conditions contain NaNs!");
            }
        }
    }

    const auto& ba      = parent->boxArray(level);
    const auto& dm      = parent->DistributionMap(level);
    const auto& ba_read = mf_read.boxArray();
    int      nc      = mf_read.nComp();

    //if we don't use a cic scheme for the initial conditions, 
    //we can safely set the number of ghost cells to 0
    //for multilevel ICs we can't use ghostcells
    mf.define(ba, dm, nc, nghost);

    mf.copy(mf_read,0,0,nc);

    if (! ((ba.contains(ba_read) && ba_read.contains(ba))) )
    {
	if (ParallelDescriptor::IOProcessor()){
            std::cout << "ba      :" << ba << std::endl;
	    std::cout << "ba_read :" << ba_read << std::endl;
            std::cout << "Read mf and hydro mf DO NOT cover the same domain!" 
        	      << std::endl;
	}
	ParallelDescriptor::Barrier();
	if (ParallelDescriptor::IOProcessor()){
            amrex::Abort();
	}
    }

    mf_read.clear();

    mf.FillBoundary();
    mf.EnforcePeriodicity(geom.periodicity());

    //FIXME
    //mf.setVal(0);

    if (mf.contains_nan())
    {
        for (int i = 0; i < mf.nComp(); i++)
        {
            if (mf.contains_nan(i, 1, nghost))
            {
                std::cout << "Found NaNs in component " << i << ". " << std::endl;
                amrex::Abort("Nyx::init_particles: Your initial conditions contain NaNs!");
            }
        }
    }
}
#endif

void Nyx::initcosmo()
{

//     if(parent->useFixedUpToLevel()<level)
//       return;
    if (ParallelDescriptor::IOProcessor())
       std::cout << "Calling InitCosmo for level " << level << std::endl;

#ifdef GRAVITY
    Real comoving_OmL;
    Real Gconst;
    const Real len[BL_SPACEDIM] = {geom.ProbLength(0),geom.ProbLength(1),geom.ProbLength(2)};

    Real particleMass;
    std::string mfDirName;

    Real redshift=-1;
    Vector<int> n_part(BL_SPACEDIM);

#ifndef NO_HYDRO
    if (level > parent->useFixedUpToLevel())
    {
        std::cout << "You have more refinement than grids, there might be a problem with your refinement criterion..." << std::endl;

    	MultiFab& S_new = get_level(level).get_new_data(State_Type);
    	MultiFab& D_new = get_level(level).get_new_data(DiagEOS_Type);

        FillCoarsePatch(S_new, 0, 0,   State_Type, 0, S_new.nComp());
        FillCoarsePatch(D_new, 0, 0, DiagEOS_Type, 0, D_new.nComp());
	return;
    }
#endif
    //
    // Read the init directory name and particle mass from the inputs file
    //
    ParmParse pp("cosmo");
    pp.get("initDirName", mfDirName);
#ifdef NUFLUID
    std::string nuMfDirName;
    pp.get("nuInitDirName", nuMfDirName);
#endif
    ParmParse pp2("nyx");
    pp2.get("initial_z",redshift);
    pp2.getarr("n_particles",n_part,0,BL_SPACEDIM);
    
#ifdef NUFLUID
    Real comoving_OmNu;
    fort_get_omnu(&comoving_OmNu);
#endif
    // fort_get_omm(&comoving_OmM );
    // fort_get_omb(&comoving_OmB );
    // fort_get_hubble(&comoving_h);
    fort_get_grav_const(&Gconst);

    // We now define this here instead of reading it.
    comoving_OmL = 1. - comoving_OmM;

    // //For whatever reason OmB is divided by OmM, which I need to undo here...
    // Real realOmB = comoving_OmB*comoving_OmM;

    //compute \Omega_K - just for checking flatness...
    Real comoving_OmK = 1 - comoving_OmM - comoving_OmL;
    if ( std::abs(comoving_OmK) > 1e-3 )
    {
            if (ParallelDescriptor::IOProcessor())
            {
        	    std::cout << "You assume a non-flat universe - \\Omega_K = "
        		      << comoving_OmK << std::endl;
            }
            amrex::Abort();
    }

    //compute \rho_{baryon}=3H_0^2\Omega_{baryon}/(8\pi G)
    Real rhoB = 3*comoving_h*100*comoving_h*100*comoving_OmB 
                / (8*M_PI*Gconst);
    if (ParallelDescriptor::IOProcessor())
    {
       std::cout << "Mean baryonic matter density is " << rhoB 
      	         << " M_sun/Mpc^3." << '\n';
    }

#ifdef NUFLUID
    //compute \rho_{\nu}=3H_0^2\Omega_{\nu}/(8\pi G)
    Real rhoNu = 3*comoving_h*100*comoving_h*100*comoving_OmNu
                 / (8*M_PI*Gconst);
    if(ParallelDescriptor::IOProcessor()){
       std::cout << "Mean neutrino matter density is " << rhoNu
                 << " M_sun/Mpc^3." << std::endl;
    }
#endif

    //compute \rho_{dark matter}=3H_0^2\Omega_{dark matter}/(8\pi G)
#ifdef NUFLUID
    Real comoving_OmD = comoving_OmM - comoving_OmB - comoving_OmNu;
#else
    Real comoving_OmD = comoving_OmM - comoving_OmB;
#endif
    Real rhoD = 3*comoving_h*100*comoving_h*100*comoving_OmD 
                / (8*M_PI*Gconst);
    if (ParallelDescriptor::IOProcessor())
    {
       std::cout << "Mean dark matter density is " << rhoD 
                  << " M_sun/Mpc^3." << '\n';
    }

    if (!do_hydro && comoving_OmB>0){
       if (ParallelDescriptor::IOProcessor()){
          std::cout << std::endl;
          std::cout << std::endl;
          std::cout << "You chose a baryonic matter content greater than 0, "
        	    << "but chose not to do hydro" << std::endl;
          std::cout << "I will be using \\rho_M for the dark matter density..." << std::endl;
          std::cout << std::endl;
          std::cout << std::endl;
       }    
       rhoD += rhoB;
    } 

    //Reads the mf and checks the data...
    MultiFab mf;
    icReadAndPrepareFab(mfDirName, 0, mf);
#ifdef NUFLUID
    MultiFab nuMf;
    icReadAndPrepareFab(nuMfDirName, 0, nuMf);
#endif

    // we have to calculate the initial a on our own
    // as there is no code path to get_comoving_a 
    // (Castro.H sources Particles.H)
    Real comoving_a = 1/(1+redshift);

    
    std::string icSource;
    pp.get("ic-source", icSource);
    int baryon_den, baryon_vx;
    int part_dx, part_vx;
    Real vel_fac[BL_SPACEDIM], dis_fac[BL_SPACEDIM];
    const Real* dx = geom.CellSize();

#ifdef NUFLUID
    int nu_den, nu_vx, nu_vy, nu_vz;
#endif

    if (icSource == "MUSIC")
    {
       if (do_hydro)
       {
          baryon_den = 0;
          baryon_vx = 1;
          part_dx = 4;
          part_vx = 7;
#ifdef NUFLUID
          nu_den = -1;
          nu_vx = -1;
#endif
       }
       else
       {
          baryon_den = -1;
          baryon_vx = -1;
          part_dx = 0;
          part_vx = 3;
#ifdef NUFLUID
          nu_den = -1;
          nu_vx = -1;
#endif
       }

#ifdef FDM
       baryon_den = 0;
       baryon_vx = 1;
       part_dx = 4;
       part_vx = 7;
#endif

       for (int n=0; n < BL_SPACEDIM; n++)
       {
          //vel_fac[n] = comoving_a * len[n]/comoving_h;
          //vel_fac[n] = len[n]/comoving_h;
          vel_fac[n] = len[n] * comoving_h; // we need the boxlength in Mpc/h
          dis_fac[n] = len[n];
       }
       particleMass = rhoD * dx[0] * dx[1] * dx[2];
       if (ParallelDescriptor::IOProcessor())
       {
          std::cout << "Particle mass at level " << level 
		    << " is " << particleMass 
                    << " M_sun." << '\n';
       }
    }
    else if (icSource == "nyx")
    {
       baryon_den = 3;
       baryon_vx = 0;
       part_dx = 0;
       part_vx = 0;
#ifdef NUFLUID
       nu_den = 3;
       nu_vx = 0;
#endif
       // Velocities are proportional to displacements by
       // LBox [MPc] * a * H [km/s/MPc]
       for (int n=0;n<BL_SPACEDIM;n++)
       {
          vel_fac[n] = len[n]*comoving_a*std::sqrt(comoving_OmM/pow(comoving_a,3)+comoving_OmL)*comoving_h*100;
	  dis_fac[n] = len[n];
       }
       //Compute particle mass
       Real simulationVolume  = len[0]*len[1]*len[2];
       Real  numberOfParticles = n_part[0] * n_part[1] * n_part[2];
       particleMass = rhoD * simulationVolume / numberOfParticles;
       if (ParallelDescriptor::IOProcessor())
       {
          std::cout << "Particle mass is " << particleMass 
                    << " M_sun." << '\n';
       }
    }
    else
    {
       std::cout << "No clue from which code the initial coniditions originate..." << std::endl
	         << "Aborting...!" << std::endl;
       amrex::Abort();
    }

    
    BoxArray baWhereNot;
    if (level < parent->initialBaLevels())
    {
//       baWhereNot = parent->boxArray(level+1);
//       baWhereNot.coarsen(parent->refRatio()[level]);
        baWhereNot = parent->initialBa(level+1);
	//std::cout << "Don't generate particles in region: " << baWhereNot << std::endl;
    }
    //initialize particles
    //std::cout << "FIXME!!!!!!!!!!!!!!!!!!" << std::endl;
    BoxArray myBaWhereNot(baWhereNot.size());
    for (int i=0; i < baWhereNot.size(); i++)
       myBaWhereNot.set(i, baWhereNot[i]);
    if (level < parent->initialBaLevels())
       myBaWhereNot.coarsen(parent->refRatio(level));
    
    // we also have to restrict the creation of particles to non refined parts of the domain.
//    if (level == 0)

    if(do_dm_particles){
      Nyx::theDMPC()->InitCosmo1ppcMultiLevel(mf, dis_fac, vel_fac, 
					      particleMass, 
					      part_dx, part_vx,
					      myBaWhereNot,
					      level, parent->initialBaLevels()+1);
    }

//    Nyx::theDMPC()->InitCosmo(mf, vel_fac, n_part, particleMass);
    //Nyx::theDMPC()->InitCosmo(mf, vel_fac, n_part, particleMass, part_dx, part_vx);

#ifdef FDM

    /*Copy density & velocity from baryon to axion state*/
    MultiFab& Ax_new = get_level(level).get_new_data(Axion_Type);
    Ax_new.setVal(0.);
    Ax_new.ParallelCopy(mf, baryon_den, 0, 4, 0, Ax_new.nGrow(),parent->Geom(level).periodicity(), FabArrayBase::COPY);
    Ax_new.plus(1, 0, 1, Ax_new.nGrow());
    Ax_new.mult(rhoD, 0, 1, Ax_new.nGrow());
    Ax_new.mult(vel_fac[0], 1, 1, Ax_new.nGrow());
    Ax_new.mult(vel_fac[1], 2, 1, Ax_new.nGrow());
    Ax_new.mult(vel_fac[2], 3, 1, Ax_new.nGrow());

    /*Once velocity is defined on all levels, construct phase from velocity divergence*/
    if(level==(parent->initialBaLevels())){                                                                          

      Vector<std::unique_ptr<MultiFab> > div(parent->initialBaLevels()+1);
      Vector<MultiFab*> phase(parent->initialBaLevels()+1);
      Vector<Vector<std::unique_ptr<MultiFab> > > grad_phase(parent->initialBaLevels()+1);
      Vector<std::array<MultiFab*,AMREX_SPACEDIM> > grad_phase_aa;

      for(int lev=0;lev<=parent->initialBaLevels();lev++){

	const BoxArray& ba = parent->boxArray()[lev];
	const DistributionMapping& dmap = parent->DistributionMap()[lev]; 
	div[lev].reset(new MultiFab(ba, dmap, 1, 0));
	div[lev]->setVal(0.0);
	phase[lev] = new MultiFab(ba, dmap, 1, 1);
	phase[lev]->setVal(0.0);
	
	MultiFab& Axvel = get_level(lev).get_new_data(Axion_Type);
	for (MFIter mfi(Axvel); mfi.isValid(); ++mfi){

	   const Box& bx = mfi.validbox();
	   const Dim3 lo = amrex::lbound(bx);
	   const Dim3 hi = amrex::ubound(bx);
	   Array4<Real> const& divarr = (*(div[lev]))[mfi].array();
	   Array4<Real> const& velarr = Axvel.array(mfi);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	   for (int k = lo.z; k <= hi.z; ++k) {
	     for (int j = lo.y; j <= hi.y; ++j) {
	       for (int i = lo.x; i <= hi.x; ++i) {
		 divarr(i,j,k,0) = (velarr(i+1,j,k,1)-velarr(i-1,j,k,1))/dx[0]*comoving_a/2.0
		   +(velarr(i,j+1,k,2)-velarr(i,j-1,k,2))/dx[1]*comoving_a/2.0
		   +(velarr(i,j,k+1,3)-velarr(i,j,k-1,3))/dx[2]*comoving_a/2.0;
	       }
	     }
	   }

	}
     
	for (int i = 0; i < (*(div[lev])).nComp(); i++)
	  {
	    if ((*(div[lev])).contains_nan(i,1,0))
	      {
		std::cout << "Testing component i for NaNs: " << i << std::endl;
		amrex::Abort("div has NaNs in this component::initcosmo");
	      }
	    if ((*(div[lev])).contains_nan(i, 1, (*(div[lev])).nGrow()))
	      {
		std::cout << "Testing component i for NaNs: " << i << std::endl;
		amrex::Abort("div has NaNs _in ghostzones_ in this component::initcosmo");
	      }
	  }

	grad_phase[lev].resize(BL_SPACEDIM);
	for (int n = 0; n < BL_SPACEDIM; ++n)
	  grad_phase[lev][n].reset(new MultiFab(BoxArray(ba).surroundingNodes(n), dmap, 1, 1));
	grad_phase_aa.push_back({AMREX_D_DECL(GetVecOfVecOfPtrs(grad_phase)[lev][0],
					      GetVecOfVecOfPtrs(grad_phase)[lev][1],
					      GetVecOfVecOfPtrs(grad_phase)[lev][2])});
      }
      
      const MultiFab* crse_bcdata = nullptr;
      Real rel_eps = 1.e-10;
      Real abs_eps = 0.;
      gravity->solve_with_MLMG(0, parent->initialBaLevels(), phase, amrex::GetVecOfConstPtrs(div),
			       grad_phase_aa, crse_bcdata, rel_eps, abs_eps);

      /*Fill FDM particle container and axion state with initial information*/
      for(int lev=0;lev<=parent->initialBaLevels();lev++){

	BoxArray baWhereNotfdm;
	if (lev < parent->initialBaLevels())
    	  baWhereNotfdm = parent->initialBa(lev+1);
    	BoxArray myBaWhereNotfdm(baWhereNotfdm.size());
    	for (int i=0; i < baWhereNotfdm.size(); i++)
    	  myBaWhereNotfdm.set(i, baWhereNotfdm[i]);
    	if (lev < parent->initialBaLevels())
    	  myBaWhereNotfdm.coarsen(parent->refRatio(lev));

	MultiFab& Ax_new = get_level(lev).get_new_data(Axion_Type);
	
	if(partlevel){
	  if(wkb_approx)
	    Nyx::theFDMwkbPC()->InitCosmo1ppcMultiLevel(Ax_new, (*phase[lev]), gamma_fdm, particleMass, myBaWhereNotfdm, lev, parent->initialBaLevels()+1);
	  else if(phase_approx)
	    Nyx::theFDMphasePC()->InitCosmo1ppcMultiLevel(Ax_new, (*phase[lev]), gamma_fdm, particleMass, myBaWhereNotfdm, lev, parent->initialBaLevels()+1);
	  else
	    amrex::Abort("Particle based FDM cosmology only works with wkb_approx=1 or phase_approx=1!");
	}

	Ax_new.ParallelCopy((*phase[lev]), 0, Nyx::AxPhas, 1, 0, Ax_new.nGrow(),parent->Geom(lev).periodicity(), FabArrayBase::COPY);
	for (MFIter mfi(Ax_new,false); mfi.isValid(); ++mfi){
	  Array4<Real> const& arr = Ax_new[mfi].array();   
	  const Box& bx = mfi.validbox();
	  const Dim3 lo = amrex::lbound(bx);
	  const Dim3 hi = amrex::ubound(bx);
	  const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};
#ifdef _OPENMP
#pragma omp parallel for
#endif
	  for(size_t i=0; i<(size_t)w.x; i++) {
	    for(size_t j=0; j<(size_t)w.y; j++) {
	      for(size_t k=0; k<(size_t)w.z; k++) {
		arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxRe)=std::sqrt(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens))*std::cos(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxPhas)/hbaroverm);
		arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxIm)=std::sqrt(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens))*std::sin(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxPhas)/hbaroverm);
	      }}}
	}
	Ax_new.FillBoundary(geom.periodicity());
      }
    }

    /*const BoxArray& ba = parent->boxArray(level);
    const DistributionMapping& dmap = parent->DistributionMap(level); 
    MultiFab vel(ba, dmap, AMREX_SPACEDIM, 1);
    vel.ParallelCopy(mf, baryon_vx, 0, AMREX_SPACEDIM, mf.nGrow(), vel.nGrow(), 
		     parent->Geom(level).periodicity(), FabArrayBase::COPY);
    vel.mult(vel_fac[0], 0, 1, vel.nGrow());
    vel.mult(vel_fac[1], 1, 1, vel.nGrow());
    vel.mult(vel_fac[2], 2, 1, vel.nGrow());

    for (int i = 0; i < vel.nComp(); i++)
      {
	if (vel.contains_nan(i,1,0))
	  {
	    std::cout << "Testing component i for NaNs: " << i << std::endl;
	    amrex::Abort("vel has NaNs in this component::initcosmo");
	  }
        if (vel.contains_nan(i, 1, vel.nGrow()))
          {
	    std::cout << "Testing component i for NaNs: " << i << std::endl;
	    amrex::Abort("vel has NaNs _in ghostzones_ in this component::initcosmo");
          }
      }
    
      MultiFab dens(ba, dmap, 1, 0);
    MultiFab::Copy(dens, mf, baryon_den, 0, 1, dens.nGrow());
    dens.plus(1, 0, 1, dens.nGrow());
    dens.mult(rhoD, 0, 1, dens.nGrow());

    for (int i = 0; i < dens.nComp(); i++)
      {
	if (dens.contains_nan(i,1,0))
	  {
	    std::cout << "Testing component i for NaNs: " << i << std::endl;
	    amrex::Abort("dens has NaNs in this component::initcosmo");
	  }
        if (dens.contains_nan(i, 1, dens.nGrow()))
          {
	    std::cout << "Testing component i for NaNs: " << i << std::endl;
	    amrex::Abort("dens has NaNs _in ghostzones_ in this component::initcosmo");
          }
      }

      MultiFab div(ba, dmap, 1, 0);
    for (MFIter mfi(vel); mfi.isValid(); ++mfi){
      const Box& bx = mfi.validbox();
      const Dim3 lo = amrex::lbound(bx);
      const Dim3 hi = amrex::ubound(bx);
      Array4<Real> const& divarr = div.array(mfi);
      Array4<Real> const& velarr = vel.array(mfi);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int k = lo.z; k <= hi.z; ++k) {
	for (int j = lo.y; j <= hi.y; ++j) {
	  for (int i = lo.x; i <= hi.x; ++i) {
	    divarr(i,j,k,0) = (velarr(i+1,j,k,0)-velarr(i-1,j,k,0))/dx[0]*comoving_a/2.0
	      +(velarr(i,j+1,k,1)-velarr(i,j-1,k,1))/dx[1]*comoving_a/2.0
	      +(velarr(i,j,k+1,2)-velarr(i,j,k-1,2))/dx[2]*comoving_a/2.0;
	  }
	}
      }

    }

    for (int i = 0; i < div.nComp(); i++)
      {
	if (div.contains_nan(i,1,0))
	  {
	    std::cout << "Testing component i for NaNs: " << i << std::endl;
	    amrex::Abort("div has NaNs in this component::initcosmo");
	  }
        if (div.contains_nan(i, 1, div.nGrow()))
          {
	    std::cout << "Testing component i for NaNs: " << i << std::endl;
	    amrex::Abort("div has NaNs _in ghostzones_ in this component::initcosmo");
          }
      }


    // divergence(mfi.validbox(),div.array(mfi),vel.array(mfi),dx,comoving_a);
    MultiFab phase(ba, dmap, 1, 1);
    Vector<std::unique_ptr<MultiFab> > grad_phase;
    Vector<std::array<MultiFab*,AMREX_SPACEDIM> > grad_phase_aa;
    grad_phase.resize(BL_SPACEDIM);
    for (int n = 0; n < BL_SPACEDIM; ++n)
      grad_phase[n].reset(new MultiFab(BoxArray(ba).surroundingNodes(n), dmap, 1, 1));
    grad_phase_aa.push_back({AMREX_D_DECL(GetVecOfPtrs(grad_phase)[0],
					  GetVecOfPtrs(grad_phase)[1],
					  GetVecOfPtrs(grad_phase)[2])});
    const MultiFab* crse_bcdata = nullptr;
    Real rel_eps = 1.e-8;
    Real abs_eps = 0.;
    gravity->solve_with_MLMG(level, level, {&phase}, {&div},
			     grad_phase_aa, crse_bcdata, rel_eps, abs_eps);
      
    for (int i = 0; i < phase.nComp(); i++)
      {
	if (phase.contains_nan(i,1,0))
	  {
	    std::cout << "Testing component i for NaNs: " << i << std::endl;
	    amrex::Abort("phase has NaNs in this component::initcosmo");
	  }
        if (phase.contains_nan(i, 1, phase.nGrow()))
          {
	    std::cout << "Testing component i for NaNs: " << i << std::endl;
	    amrex::Abort("phase has NaNs _in ghostzones_ in this component::initcosmo");
          }
      }

    if (phase.contains_nan())
      {
	std::cout << "Found NaNs in initial FDM phase. " << std::endl;
	amrex::Abort("Nyx::init_particles: Your initial FDM phase contain NaNs!");
      }


    if(partlevel){
      if(wkb_approx)
	Nyx::theFDMwkbPC()->InitCosmo1ppcMultiLevel(vel, phase, dens, gamma_fdm, particleMass, myBaWhereNot, level, parent->initialBaLevels()+1);
      else if(phase_approx)
	Nyx::theFDMphasePC()->InitCosmo1ppcMultiLevel(vel, phase, dens, gamma_fdm, particleMass, myBaWhereNot, level, parent->initialBaLevels()+1);
      else
	amrex::Abort("FDM cosmology only works with wkb_approx=1 or phase_approx=1!");
      // Nyx::theFDMPC()->InitCosmo1ppcMultiLevel(particle_mf, phase, gamma_fdm, particleMass, myBaWhereNotfdm, lev, parent->initialBaLevels()+1);
      //}

      //    if(partlevel && level==(parent->initialBaLevels())){                                                                                                      
      if(level==(parent->initialBaLevels())){                                                                          
                          
       Vector<std::unique_ptr<MultiFab> > particle_mf;
       Vector<std::unique_ptr<MultiFab> > div(parent->initialBaLevels()+1);
       Vector<MultiFab*> phase(parent->initialBaLevels()+1);
       Vector<Vector<std::unique_ptr<MultiFab> > > grad_phase(parent->initialBaLevels()+1);
       Vector<std::array<MultiFab*,AMREX_SPACEDIM> > grad_phase_aa;

       if(wkb_approx)
	 Nyx::theFDMwkbPC()->AssignDensityAndVels(particle_mf);
       else if(phase_approx)
	 Nyx::theFDMphasePC()->AssignDensityAndVels(particle_mf);
       else
	 amrex::Abort("FDM cosmology only works with wkb_approx=1 or phase_approx=1!");
       

       for(int lev=0;lev<=parent->initialBaLevels();lev++){

	 const BoxArray& ba = parent->boxArray()[lev];
	 const DistributionMapping& dmap = parent->DistributionMap()[lev]; 
	 div[lev].reset(new MultiFab(ba, dmap, 1, 0));
	 div[lev]->setVal(0.0);
	 phase[lev] = new MultiFab(ba, dmap, 1, 1);
	 phase[lev]->setVal(0.0);

	 for (MFIter mfi(*(particle_mf[lev])); mfi.isValid(); ++mfi){
	   FArrayBox&  vel  = (*(particle_mf[lev]))[mfi];
	   const Box& bx = mfi.validbox();
	   const Dim3 lo = amrex::lbound(bx);
	   const Dim3 hi = amrex::ubound(bx);
	   Array4<Real> const& divarr = div.array(mfi);
	   Array4<Real> const& velarr = vel.array(mfi);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	   for (int k = lo.z; k <= hi.z; ++k) {
	     for (int j = lo.y; j <= hi.y; ++j) {
	       for (int i = lo.x; i <= hi.x; ++i) {
		 divarr(i,j,k,0) = (velarr(i+1,j,k,1)-velarr(i-1,j,k,1))/dx[0]*comoving_a/2.0
		   +(velarr(i,j+1,k,2)-velarr(i,j-1,k,2))/dx[1]*comoving_a/2.0
		   +(velarr(i,j,k+1,3)-velarr(i,j,k-1,3))/dx[2]*comoving_a/2.0;
	       }
	     }
	   }

	 }

	 for (int i = 0; i < div.nComp(); i++)
	   {
	     if (div.contains_nan(i,1,0))
	       {
		 std::cout << "Testing component i for NaNs: " << i << std::endl;
		 amrex::Abort("div has NaNs in this component::initcosmo");
	       }
	     if (div.contains_nan(i, 1, div.nGrow()))
	       {
		 std::cout << "Testing component i for NaNs: " << i << std::endl;
		 amrex::Abort("div has NaNs _in ghostzones_ in this component::initcosmo");
	       }
	   }
	 



	 //for (MFIter mfi(*(particle_mf[lev])); mfi.isValid(); ++mfi){
	 //  FArrayBox&  vel  = (*(particle_mf[lev]))[mfi];
	 //  const Box& box = mfi.tilebox();
	 //  const int* lo = box.loVect();
	 //  const int* hi = box.hiVect();
	 //  const Real* dx_lev = parent->Geom(lev).CellSize();
	 //  Array4<Real> const& div_array = (*(div[lev]))[mfi].array();
	 //  Array4<Real const> const& vel_array = vel.array();
	 //  divergence(box,div_array,vel_array,dx_lev,comoving_a);
	 //}
	 
	 grad_phase[lev].resize(BL_SPACEDIM);
	 for (int n = 0; n < BL_SPACEDIM; ++n)
	   grad_phase[lev][n].reset(new MultiFab(BoxArray(ba).surroundingNodes(n), dmap, 1, 1));
	 grad_phase_aa.push_back({AMREX_D_DECL(GetVecOfVecOfPtrs(grad_phase)[lev][0],
					       GetVecOfVecOfPtrs(grad_phase)[lev][1],
					       GetVecOfVecOfPtrs(grad_phase)[lev][2])});
       }

       const MultiFab* crse_bcdata = nullptr;
       Real rel_eps = 1.e-10;
       Real abs_eps = 0.;
       gravity->solve_with_MLMG(0, parent->initialBaLevels(), phase, amrex::GetVecOfConstPtrs(div),
				grad_phase_aa, crse_bcdata, rel_eps, abs_eps);

       for(int lev=0;lev<=parent->initialBaLevels();lev++){
	 BoxArray baWhereNotfdm;
	 if (lev < parent->initialBaLevels())
    //	  baWhereNotfdm = parent->initialBa(lev+1);
    //	BoxArray myBaWhereNotfdm(baWhereNotfdm.size());
    //	for (int i=0; i < baWhereNotfdm.size(); i++)
    //	  myBaWhereNotfdm.set(i, baWhereNotfdm[i]);
    //	if (lev < parent->initialBaLevels())
    //	  myBaWhereNotfdm.coarsen(parent->refRatio(lev));

    //if(wkb_approx)
    //Nyx::theFDMwkbPC()->InitCosmo1ppcMultiLevel(particle_mf, phase, gamma_fdm, particleMass, myBaWhereNotfdm, lev, parent->initialBaLevels()+1);
    //else if(phase_approx)
    //Nyx::theFDMphasePC()->InitCosmo1ppcMultiLevel(particle_mf, phase, gamma_fdm, particleMass, myBaWhereNotfdm, lev, parent->initialBaLevels()+1);
    //else
    //amrex::Abort("FDM cosmology only works with wkb_approx=1 or phase_approx=1!");
	// Nyx::theFDMPC()->InitCosmo1ppcMultiLevel(particle_mf, phase, gamma_fdm, particleMass, myBaWhereNotfdm, lev, parent->initialBaLevels()+1);
    //}
    //}

//for(int lev=0;lev<=parent->initialBaLevels();lev++){
    MultiFab& Ax_new = get_level(level).get_new_data(Axion_Type);
    Ax_new.setVal(0.);

    Ax_new.ParallelCopy(dens, 0, Nyx::AxDens, 1, 0, Ax_new.nGrow(),parent->Geom(level).periodicity(), FabArrayBase::COPY);
    Ax_new.ParallelCopy(phase, 0, Nyx::AxPhas, 1, 0, Ax_new.nGrow(),parent->Geom(level).periodicity(), FabArrayBase::COPY);

    for (MFIter mfi(Ax_new,false); mfi.isValid(); ++mfi){
      Array4<Real> const& arr = Ax_new[mfi].array();
      //      Array4<Real> const& axold = Ax_old[mfi].array();                                                                                                           
      const Box& bx = mfi.validbox();
      const Dim3 lo = amrex::lbound(bx);
      const Dim3 hi = amrex::ubound(bx);
      const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(size_t i=0; i<(size_t)w.x; i++) {
	for(size_t j=0; j<(size_t)w.y; j++) {
	  for(size_t k=0; k<(size_t)w.z; k++) {
	    //              size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*i+(size_t)w.z*j+k;                                                                      
	    //complex_t temp = a[local_indx_threaded];                                                                                                                   
	    arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxRe)=std::sqrt(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens))*std::cos(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxPhas)/hbaroverm);
	    arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxIm)=std::sqrt(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens))*std::sin(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxPhas)/hbaroverm);
	    //              arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens)=std::real(temp)*std::real(temp)+std::imag(temp)*std::imag(temp);                                       
	    //arr(i+lo.x,j+lo.y,k+lo.z,Nyx::AxPhas) = std::arg(temp);                                                                                                    
	  }}}
    }//}

//MultiFab& Ax_new = get_new_data(Axion_Type);
    Ax_new.FillBoundary(geom.periodicity());

    //}

    //MultiFab& Ax_new = get_new_data(Axion_Type);
    //Ax_new.setVal(0.);


    // if(partlevel){
      
    //   const BoxArray& ba = parent->boxArray(level);
    //   const DistributionMapping& dmap = parent->DistributionMap(level); 
    //   MultiFab vel(ba, dmap, AMREX_SPACEDIM, 1);
    //   vel.ParallelCopy(mf, baryon_vx, 0, AMREX_SPACEDIM, mf.nGrow(), vel.nGrow(), 
    // 		       parent->Geom(level).periodicity(), FabArrayBase::COPY);
    //   vel.mult(vel_fac[0], 0, 1, vel.nGrow());
    //   vel.mult(vel_fac[1], 1, 1, vel.nGrow());
    //   vel.mult(vel_fac[2], 2, 1, vel.nGrow());

    //   MultiFab dens(ba, dmap, 1, 0);
    //   MultiFab::Copy(dens, mf, baryon_den, 0, 1, dens.nGrow());
    //   dens.plus(1, 0, 1, dens.nGrow());
    //   dens.mult(rhoD, 0, 1, dens.nGrow());
      
    //   MultiFab div(ba, dmap, 1, 0);
    //   for (MFIter mfi(vel); mfi.isValid(); ++mfi)
    // 	divergence(mfi.validbox(),div.array(mfi),vel.array(mfi),dx,comoving_a);
    //   MultiFab phase(ba, dmap, 1, 1);
    //   Vector<std::unique_ptr<MultiFab> > grad_phase;
    //   Vector<std::array<MultiFab*,AMREX_SPACEDIM> > grad_phase_aa;
    //   grad_phase.resize(BL_SPACEDIM);
    //   for (int n = 0; n < BL_SPACEDIM; ++n)
    // 	grad_phase[n].reset(new MultiFab(BoxArray(ba).surroundingNodes(n), dmap, 1, 1));
    //   grad_phase_aa.push_back({AMREX_D_DECL(GetVecOfPtrs(grad_phase)[0],
    // 					    GetVecOfPtrs(grad_phase)[1],
    // 					    GetVecOfPtrs(grad_phase)[2])});
    //   const MultiFab* crse_bcdata = nullptr;
    //   Real rel_eps = 1.e-8;
    //   Real abs_eps = 0.;
    //   gravity->solve_with_MLMG(level, level, {&phase}, {&div},
    // 		      grad_phase_aa, crse_bcdata, rel_eps, abs_eps);
      
    //   if (phase.contains_nan())
    // 	{
    // 	  std::cout << "Found NaNs in initial FDM phase. " << std::endl;
    // 	  amrex::Abort("Nyx::init_particles: Your initial FDM phase contain NaNs!");
    // 	}
      
    //   BoxArray baWhereNotfdm;
    //   if (level < parent->initialBaLevels())
    // 	baWhereNotfdm = parent->initialBa(level+1);
    //   BoxArray myBaWhereNotfdm(baWhereNotfdm.size());
    //   for (int i=0; i < baWhereNotfdm.size(); i++)
    // 	myBaWhereNotfdm.set(i, baWhereNotfdm[i]);
    //   if (level < parent->initialBaLevels())
    // 	myBaWhereNotfdm.coarsen(parent->refRatio(level));
    //   if(wkb_approx)
    // 	Nyx::theFDMwkbPC()->InitCosmo1ppcMultiLevel(vel, phase, dens, gamma_fdm, particleMass, myBaWhereNotfdm, 
    // 						    level, parent->initialBaLevels()+1);
    //   else
    // 	amrex::Abort("FDM cosmology only works with wkb_approx=1!");
    // }

    // MultiFab& Ax_new = get_new_data(Axion_Type);
    // Ax_new.setVal(0.);
    */
#endif

#ifndef NO_HYDRO
    MultiFab& S_new = get_level(level).get_new_data(State_Type);
    MultiFab& D_new = get_level(level).get_new_data(DiagEOS_Type);

#ifdef NUFLUID
    //copy density 
    S_new.copy(nuMf, nu_den, NuDens, 1);
    S_new.plus(1,            NuDens, 1, S_new.nGrow());
    S_new.mult(rhoNu,        NuDens, 1, S_new.nGrow());

    //copy velocities...
    S_new.copy(nuMf, nu_vx,  NuXMom, 3);

    //...and "transform" to momentum
    S_new.mult(vel_fac[0],   NuXMom, 1, S_new.nGrow());
    S_new.mult(vel_fac[1],   NuYMom, 1, S_new.nGrow());
    S_new.mult(vel_fac[2],   NuZMom, 1, S_new.nGrow());
    MultiFab::Multiply(S_new, S_new, NuDens, NuXMom, 1, S_new.nGrow());
    MultiFab::Multiply(S_new, S_new, NuDens, NuYMom, 1, S_new.nGrow());
    MultiFab::Multiply(S_new, S_new, NuDens, NuZMom, 1, S_new.nGrow());
#endif


    //initialize hydro
    //units for velocity should be okay
    if (do_hydro)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "Do hydro initialization..." << '\n';
        }

    	MultiFab& S_new = get_level(level).get_new_data(State_Type);
    	MultiFab& D_new = get_level(level).get_new_data(DiagEOS_Type);

        //Fill everything with old data...should only affect ghostzones, but
	//seems to have no effect...
	if (level > 0)
	{
           FillCoarsePatch(S_new, 0, 0,   State_Type, 0, S_new.nComp());
	   FillCoarsePatch(D_new, 0, 0, DiagEOS_Type, 0, D_new.nComp());
	}

     	//copy density 
     	S_new.setVal(0.);
     	S_new.copy(mf, baryon_den, Density, 1);
     	S_new.plus(1,     Density, 1, S_new.nGrow());
     	S_new.mult(rhoB,  Density, 1, S_new.nGrow());

//      //This block assigns "the same" density for the baryons as for the dm.
//      Vector<std::unique_ptr<MultiFab> > particle_mf;
//      Nyx::theDMPC()->AssignDensity(particle_mf);
//      particle_mf[0]->mult(comoving_OmB / comoving_OmD);
//      S_new.copy(*particle_mf[0], 0, Density, 1);

     	//copy velocities...
     	S_new.copy(mf, baryon_vx, Xmom, 3);

        //...and "transform" to momentum
     	S_new.mult(vel_fac[0], Xmom, 1, S_new.nGrow());
     	MultiFab::Multiply(S_new, S_new, Density, Xmom, 1, S_new.nGrow());
     	S_new.mult(vel_fac[1], Ymom, 1, S_new.nGrow());
     	MultiFab::Multiply(S_new, S_new, Density, Ymom, 1, S_new.nGrow());
     	S_new.mult(vel_fac[2], Zmom, 1, S_new.nGrow());
        MultiFab::Multiply(S_new, S_new, Density, Zmom, 1, S_new.nGrow());

        // Convert (rho X)_i to X_i before calling init_e_from_T
//      if (use_const_species == 0)
//          for (int i = 0; i < NumSpec; i++) 
//              MultiFab::Divide(S_new, S_new, Density, FirstSpec+i, 1, 0);

        Real tempInit = 0.021*(1.0+redshift)*(1.0+redshift);

        int ns = S_new.nComp();
        int nd = D_new.nComp();

        D_new.setVal(tempInit, Temp_comp);
        D_new.setVal(0.0, Ne_comp);
        if (inhomo_reion > 0)
            D_new.setVal(0.0, Zhi_comp);

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
      	        BL_TO_FORTRAN(D_new[mfi]), &nd, lo, hi, &old_a);
     	}     	

        // Convert X_i to (rho X)_i
        if (use_const_species == 0)
        {
           S_new.setVal(0.75, FirstSpec);
           S_new.setVal(0.25, FirstSpec+1);

           for (int i = 0; i < NumSpec; i++)
           {
              MultiFab::Multiply(S_new, S_new, Density, FirstSpec+i, 1, 0);
           }
        }

    }
    else
    {
       // S_new.setVal(0.0, Density);
       // D_new.setVal(-23, Temp_comp);
       // D_new.setVal(-42, Ne_comp);
    }
#endif //ifndef NO_HYDRO       
    mf.clear();

#endif //ifdef GRAVITY
}
