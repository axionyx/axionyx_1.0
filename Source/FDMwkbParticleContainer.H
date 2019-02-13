#ifdef FDM
#ifndef _FDMwkbParticleContainer_H_
#define _FDMwkbParticleContainer_H_

#include "NyxParticleContainer.H"

class FDMwkbParticleContainer
    : public NyxParticleContainer<26>
{
public:
    FDMwkbParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<26>(amr)
    {
      real_comp_names.clear();
      real_comp_names.push_back("mass");
      real_comp_names.push_back("xvel");
      real_comp_names.push_back("yvel");
      real_comp_names.push_back("zvel");
      real_comp_names.push_back("phase");
      real_comp_names.push_back("amp1");
      real_comp_names.push_back("amp2");
      real_comp_names.push_back("width");
      real_comp_names.push_back("qq1");
      real_comp_names.push_back("qq2");
      real_comp_names.push_back("qq3");
      real_comp_names.push_back("qq4");
      real_comp_names.push_back("qq5");
      real_comp_names.push_back("qq6");
      real_comp_names.push_back("qq7");
      real_comp_names.push_back("qq8");
      real_comp_names.push_back("qq9");
      real_comp_names.push_back("pq1");
      real_comp_names.push_back("pq2");
      real_comp_names.push_back("pq3");
      real_comp_names.push_back("pq4");
      real_comp_names.push_back("pq5");
      real_comp_names.push_back("pq6");
      real_comp_names.push_back("pq7");
      real_comp_names.push_back("pq8");
      real_comp_names.push_back("pq9");
    }

    using MyParIter = amrex::ParIter<26>;
    using MyConstParIter = amrex::ParConstIter<26>;

    virtual ~FDMwkbParticleContainer () {}

    void InitCosmo  (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Vector<int> n_part, const amrex::Real particleMass);
    void InitCosmo  (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Vector<int> n_part,
                     const amrex::Real particleMass, const amrex::Real shift[]);
    void InitCosmo1ppc (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Real particleMass);
    void InitCosmo1ppcMultiLevel(amrex::MultiFab& mf, const amrex::Real disp_fac[], const amrex::Real vel_fac[], 
                                 const amrex::Real particleMass, int disp_idx, int vel_idx, 
                                 amrex::BoxArray &baWhereNot, int lev, int nlevs);

    void CreateGhostParticlesFDM (int level, int lev, int nGrow, AoS& ghosts) const;

    void DepositFDMParticles (amrex::MultiFab& mf_real, amrex::MultiFab& mf_imag, int level, amrex::Real a) const;

    amrex::Real estTimestepFDM(amrex::MultiFab& phi, int lev, amrex::Real cfl) const;

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

    void InitGaussianBeams (long n_axpart, int lev, int nlevs, const amrex::Real hbaroverm, const amrex::Real sigma_ax, const amrex::Real gamma_ax);

    amrex::Real generateGaussianNoise(const amrex::Real &mean, const amrex::Real &stdDev);
};

#endif /* _FDMwkbParticleContainer_H_ */
#endif /* FDM */
