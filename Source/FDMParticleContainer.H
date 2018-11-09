#ifdef FDM
#ifndef _FDMParticleContainer_H_
#define _FDMParticleContainer_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include "NyxParticleContainer.H"

// We make FDMParticleContainer a class instead of a typedef so that
//    we can have Nyx-specific functions here instead of in amrex
class FDMParticleContainer
  : public NyxParticleContainer<44, 0, 0, 0>
{

    using MyParIter = amrex::ParIter<44>;
    using MyConstParIter = amrex::ParConstIter<44>;

public:
    FDMParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<44>(amr)
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
      real_comp_names.push_back("qp1");
      real_comp_names.push_back("qp2");
      real_comp_names.push_back("qp3");
      real_comp_names.push_back("qp4");
      real_comp_names.push_back("qp5");
      real_comp_names.push_back("qp6");
      real_comp_names.push_back("qp7");
      real_comp_names.push_back("qp8");
      real_comp_names.push_back("qp9");
      real_comp_names.push_back("pq1");
      real_comp_names.push_back("pq2");
      real_comp_names.push_back("pq3");
      real_comp_names.push_back("pq4");
      real_comp_names.push_back("pq5");
      real_comp_names.push_back("pq6");
      real_comp_names.push_back("pq7");
      real_comp_names.push_back("pq8");
      real_comp_names.push_back("pq9");
      real_comp_names.push_back("pp1");
      real_comp_names.push_back("pp2");
      real_comp_names.push_back("pp3");
      real_comp_names.push_back("pp4");
      real_comp_names.push_back("pp5");
      real_comp_names.push_back("pp6");
      real_comp_names.push_back("pp7");
      real_comp_names.push_back("pp8");
      real_comp_names.push_back("pp9");
    }

    virtual ~FDMParticleContainer () {}

    void InitCosmo  (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Vector<int> n_part, const amrex::Real particleMass);
    void InitCosmo  (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Vector<int> n_part,
                     const amrex::Real particleMass, const amrex::Real shift[]);
    void InitCosmo1ppc (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Real particleMass);
    void InitCosmo1ppcMultiLevel(amrex::MultiFab& mf, const amrex::Real disp_fac[], const amrex::Real vel_fac[], 
                                 const amrex::Real particleMass, int disp_idx, int vel_idx, 
                                 amrex::BoxArray &baWhereNot, int lev, int nlevs);

    void AssignDensityAndVels (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0) const;

    virtual void moveKickDrift (amrex::MultiFab& acceleration, int level, amrex::Real timestep,
                                amrex::Real a_old = 1.0, amrex::Real a_half = 1.0, int where_width = 0);
    virtual void moveKick      (amrex::MultiFab& acceleration, int level, amrex::Real timestep,
                                amrex::Real a_new = 1.0, amrex::Real a_half = 1.0);

    void InitFromBinaryMortonFile(const std::string& particle_directory, int nextra, int skip_factor);
 
    void InitVarCount (amrex::MultiFab& mf, long n_axpart);
};

#endif /* _FDMParticleContainer_H_ */
#endif /* FDM */
