#ifdef FDM
#ifndef _FDMphaseParticleContainer_H_
#define _FDMphaseParticleContainer_H_

#include "NyxParticleContainer.H"

class FDMphaseParticleContainer
    : public NyxParticleContainer<5>
{
public:
    FDMphaseParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<5>(amr)
    {
      real_comp_names.clear();
      real_comp_names.push_back("mass");
      real_comp_names.push_back("xvel");
      real_comp_names.push_back("yvel");
      real_comp_names.push_back("zvel");
      real_comp_names.push_back("phase");
    }

    using MyParIter = amrex::ParIter<5>;
    using MyConstParIter = amrex::ParConstIter<5>;

    virtual ~FDMphaseParticleContainer () {}

    void InitCosmo  (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Vector<int> n_part, const amrex::Real particleMass);
    void InitCosmo  (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Vector<int> n_part,
                     const amrex::Real particleMass, const amrex::Real shift[]);
    void InitCosmo1ppc (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Real particleMass);
  void InitCosmo1ppcMultiLevel(amrex::Vector<std::unique_ptr<amrex::MultiFab> >&  mf,  amrex::Vector<amrex::MultiFab*>& phase,
  			       const amrex::Real gamma_fdm, const amrex::Real particleMass, 
  			       amrex::BoxArray &baWhereNot, int lev, int nlevs);

  void InitCosmo1ppcMultiLevel(amrex::MultiFab& vel, amrex::MultiFab& phase, amrex::MultiFab& dens,
			       const amrex::Real gamma_fdm, const amrex::Real particleMass, 
			       amrex::BoxArray &baWhereNot, int lev, int nlevs);

  void InitCosmo1ppcMultiLevel(amrex::MultiFab& densvel, amrex::MultiFab& phase,
			       const amrex::Real gamma_fdm, const amrex::Real particleMass, 
			       amrex::BoxArray &baWhereNot, int lev, int nlevs);

    void CreateGhostParticlesFDM (int level, int lev, int nGrow, AoS& ghosts) const;

  void DepositFDMParticles (amrex::MultiFab& mf_real, amrex::MultiFab& mf_imag, int level, amrex::Real a, amrex::Real theta_fdm, amrex::Real hbaroverm) const;

  void DepositFDMParticlesCWA (amrex::MultiFab& mf_real, amrex::MultiFab& mf_imag, int level, amrex::Real a, amrex::Real theta_fdm, amrex::Real hbaroverm) const;

    amrex::Real estTimestepFDM(amrex::MultiFab& phi, amrex::Real a, int lev, amrex::Real cfl) const;

    virtual void moveKickDrift (amrex::MultiFab& acceleration, int level, amrex::Real timestep,
                                amrex::Real a_old = 1.0, amrex::Real a_half = 1.0, int where_width = 0);
    virtual void moveKick      (amrex::MultiFab& acceleration, int level, amrex::Real timestep,
                                amrex::Real a_new = 1.0, amrex::Real a_half = 1.0);

    void moveKickDriftFDM (amrex::MultiFab& phi, int grav_n_grow, amrex::MultiFab& acceleration, int level, amrex::Real timestep,
			   amrex::Real a_old = 1.0, amrex::Real a_half = 1.0, int where_width = 0);
    void moveKickFDM      (amrex::MultiFab& phi, int grav_n_grow, amrex::MultiFab& acceleration, int level, amrex::Real timestep,
			   amrex::Real a_new = 1.0, amrex::Real a_half = 1.0);

    void InitFromBinaryMortonFile(const std::string& particle_directory, int nextra, int skip_factor);
 
    void InitVarCount (amrex::MultiFab& mf, long n_axpart, amrex::BoxArray &baWhereNot, int lev, int nlevs);

    void InitCWA (long n_axpart, int lev, int nlevs, const amrex::Real hbaroverm, const amrex::Real a);

};

#endif /* _FDMphaseParticleContainer_H_ */
#endif /* FDM */
