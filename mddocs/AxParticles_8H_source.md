
# File AxParticles.H

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**AxParticles.H**](AxParticles_8H.md)

[Go to the documentation of this file.](AxParticles_8H.md) 


````cpp
#ifndef _AXPARTICLES_H_
#define _AXPARTICLES_H_ 

#include <map>
#include <deque>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>

#include <string>
#include <sstream>

#include <ParmParse.H>

#include <REAL.H>
#include <IntVect.H>
#include <Array.H>
#include <Amr.H>
#include <AmrLevel.H>
#include <Utility.H>
#include <Geometry.H>
#include <VisMF.H>
#include <Particles_F.H>
#include <RealBox.H>

#include <Particles.H>

#include <AxParticles_F.H>

#ifdef _OPENMP
#include <omp.h>
#endif


template <int N>
class AxionParticleContainer
    :
    public ParticleContainer<N>
{
public:

    AxionParticleContainer (Amr* amr)
       :
       ParticleContainer<N>(amr) {}

    //copied typedefs, dont know another way to do it.
    //
    // The type of Particles we hold.
    //
    typedef Particle<N> ParticleType;
    //
    // We want to store the particles on a level by level and grid by grid basis.  This will
    // make accessing them and doing operations on them more memory efficient since most of our
    // operations on particles are done on a level by level basis or grid by grid basis.
    //
    typedef typename std::deque<ParticleType> PBox;
    //
    // A level of particles is stored in a map indexed by the grid number.
    //
    typedef typename std::map<int,PBox> PMap;

    void InitEmpty(const int max_level);

    void moveKickDrift (MultiFab& grav_vector, int level, Real timestep, 
                        Real a_old = 1.0, Real a_half = 1.0);
    void moveKick      (MultiFab& grav_vector, int level, Real timestep, 
                        Real a_new = 1.0, Real a_half = 1.0,
                        int start_comp_for_accel = -1); 
    
    Real estTimestep (MultiFab& grav_vector, Real a, int lev, Real cfl);
    
    void AssignDensity (PArray<MultiFab>& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1, double masstreshold=1.e50) const;
    void AssignDensitySingleLevel (MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const;
    
    void AssignDensityGradientsSingleLevel (MultiFab& mf_to_be_filled,
                                                int       lev) const;
    void quantummoveKick(Real timestep,Real a_old,Real a_half, bool before_drift);

    void doFiltering(MultiFab& mf,int comp, int niter, int nstride) const;

    void averageVelocity (int level);
    void remap();

    void addQuantumPressure (MultiFab & acc_vector, int lev, Real a_old, Real a_half);
    void computeQuantumPressure (MultiFab & axPressForceMf, MultiFab & pressureMf, int lev);
    void addArtificialViscosity (MultiFab & acc_vector, int lev, Real a_old, Real a_half);

    void computeEnergies(Real cur_time, MultiFab& Phi, Real a);

    void InitVarMass (MultiFab& mf);
    void InitVarCount (MultiFab& mf, long n_axpart);
    static Real max_acceleration;
    static int smoothing_length;
    static Real damping_constant;
    
    static bool particles_moved;
    static bool initzd;
    
    static MultiFab* mf_pointer;

protected:

};

template<int N> int AxionParticleContainer<N>::smoothing_length = 4;
template<int N> Real AxionParticleContainer<N>::max_acceleration = 1e20;
template<int N> Real AxionParticleContainer<N>::damping_constant = 0.0;
template<int N> bool AxionParticleContainer<N>::particles_moved = true;
template<int N> bool AxionParticleContainer<N>::initzd = false;
template<int N> MultiFab* AxionParticleContainer<N>::mf_pointer = 0;

template <int N>
void AxionParticleContainer<N>::remap(){
    
    //Assign density
    MultiFab mf(this->m_amr->ParticleBoxArray(0), 4, 1);
    ParticleContainer<N>::AssignDensitySingleLevel(mf, 0, 4);
    
    //Delete particles
    for (int lev = 0; lev < this->m_particles.size(); lev++)
    {
        this->m_particles[lev].clear();
    }
    
    //Reinitialize them
    InitVarMass(mf);
}


template <int N>
void AxionParticleContainer<N>::computeEnergies(Real cur_time, MultiFab& Phi, Real a){
    const Geometry& geom = this->m_amr->Geom(0);
    const Real* dx       = geom.CellSize();
    
    PArray<MultiFab> PMF;
    ParticleContainer<N>::AssignDensity(PMF,0,4); //non-filtered density field
    
    
    MultiFab gradkinen(this->m_amr->ParticleBoxArray(0), 1, 0);
    
    //Add one ghost cell to the non-filtered density
    MultiFab poten(this->m_amr->ParticleBoxArray(0), 1, 1);   
    poten.copy(PMF[0],0,0,1);
    if(!geom.isAllPeriodic())
      poten.setBndry(0.0);
    poten.FillBoundary();
    geom.FillPeriodicBoundary(poten);
    


         
    //Compute the densitygradient part of the kinetic energy and store it in gradkinen
    for (MFIter mfi(poten); mfi.isValid(); ++mfi)
    {
         const Box& box = mfi.validbox();
         const int*  lo = box.loVect();
         const int*  hi = box.hiVect();

         BL_FORT_PROC_CALL(KINENERGY, kinenergy)
                           (BL_TO_FORTRAN(poten[mfi]),
               BL_TO_FORTRAN(gradkinen[mfi]),
                           dx, lo, hi);
    }

    //Multiply non-filtered density by the gravitational potential to get the potential energy
    MultiFab::Multiply(poten,Phi,0,0,1,0);
    poten.mult(0.5/a);
    
    //Compute kinetic energy stored in the average movement per cell
    MultiFab::Multiply(PMF[0],PMF[0],1,1,1,0);
    MultiFab::Multiply(PMF[0],PMF[0],2,2,1,0);
    MultiFab::Multiply(PMF[0],PMF[0],3,3,1,0);
    MultiFab::Add(PMF[0],PMF[0],2,1,1,0);
    MultiFab::Add(PMF[0],PMF[0],3,1,1,0);
    MultiFab::Multiply(PMF[0],PMF[0],1,0,1,0);
    PMF[0].mult(0.5);
 
    //Sum up particle kinetic energy and momenta
    PMap&      pmap          = this->m_particles[0];
    Real particlekinen[4]    = {0.0,0.0,0.0,0.0};
    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);
        // Uncomment to compute the quantities only in a given radius 
        //   if(sqrt((p.m_pos[0]-0.128)*(p.m_pos[0]-0.128)+(p.m_pos[1]-0.128)*(p.m_pos[1]-0.128)+(p.m_pos[2]-0.128)*(p.m_pos[2]-0.128)) < 0.025){
               particlekinen[0] += 0.5 * p.m_data[0] * (p.m_data[1]*p.m_data[1] + p.m_data[2]*p.m_data[2] + p.m_data[3]*p.m_data[3]);
               for(int j=1;j<4;j++)
               particlekinen[j] += p.m_data[0]*p.m_data[j];
        //    }
        } 
    }
    ParallelDescriptor::ReduceRealSum(particlekinen,4);
    
    //Integrate potential energy, gradkinen and avkinen
    Real totalpoten = poten.norm1(0,0)*dx[0]*dx[1]*dx[2];
    Real totalgradkinen =  gradkinen.norm1(0,0)*dx[0]*dx[1]*dx[2]/a/a;
    Real totalavkinen = PMF[0].norm1(0,0)*dx[0]*dx[1]*dx[2];
    if(ParallelDescriptor::IOProcessor()){
        std::ofstream endiagfile;
    endiagfile.open("energydiag", std::ofstream::out | std::ofstream::app);
    endiagfile.precision(15);
    endiagfile << cur_time << ' ' << totalpoten << ' ' << totalgradkinen << ' ' << particlekinen[0] << " " << totalavkinen << " " << a << " " << particlekinen[1] << " " << particlekinen[2] << " " << particlekinen[3] <<'\n';
        endiagfile.close();
    }
}


template <int N>
void AxionParticleContainer<N>::InitEmpty(const int max_level){
    this->m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.
    this->m_particles.resize(max_level+1);
    for (int lev = 0; lev < this->m_particles.size(); lev++)
    {
        BL_ASSERT(this->m_particles[lev].empty());
    }
}

template <int N>
void AxionParticleContainer<N>::averageVelocity (int level)
{
    const Geometry& geom = this->m_amr->Geom(0);
    const Real* dx       = geom.CellSize();
    const Real* plo      = geom.ProbLo();
    
    PMap&      pmap          = this->m_particles[0];
    
    MultiFab mf(this->m_amr->ParticleBoxArray(0), 4, 1);
    
     //use the AssignDensity version without filterin
    ParticleContainer<N>::AssignDensitySingleLevel(mf, 0, 4); //use the
    
    //Multiply each velocity component with the density
    for( int i = 1; i<4; i++)
        MultiFab::Multiply(mf,mf,0,i,1,1);

    //Filter the momenta
    for( int i = 0; i<4; i++)
        this->doFiltering(mf,i,2,2);
    
    //Assign the velocity in each cell to the particles it contains
    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = mf[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);
            //
            // note: m_data[0] is mass, 1 is v_x, ...
            //
        Real dens;
            Real vel[BL_SPACEDIM];
            IntVect cell (  D_DECL(floor((p.m_pos[0]-plo[0])/dx[0]),
                                           floor((p.m_pos[1]-plo[1])/dx[1]),
                                           floor((p.m_pos[2]-plo[2])/dx[2])) );
            if (!gfab.box().contains(cell)) continue;
    
            dens = gfab(cell,0);
            for(int j = 0; j<BL_SPACEDIM; j++)
                   vel[j] = gfab(cell,j+1);
        
            D_TERM(p.m_data[1] = vel[0]/dens;,
                   p.m_data[2] = vel[1]/dens;,
                   p.m_data[3] = vel[2]/dens;);
        } 
      
    }
}
    
template <int N>
void AxionParticleContainer<N>::moveKickDrift (MultiFab& grav_vector, 
                int level, Real timestep, 
                    Real a_old, Real a_half)
{
   // grav_vector.setVal(0.0);
    quantummoveKick(timestep, a_half, a_old,true);
   
    ParticleContainer<N>::moveKickDrift(grav_vector, level, timestep, a_old, a_half);
    ParticleContainer<N>::Redistribute(false,true);
    particles_moved = true;
}

template <int N>
void AxionParticleContainer<N>::moveKick (MultiFab& grav_vector, 
            int level, Real timestep, 
                    Real a_new, Real a_half,
                    int start_comp_for_accel)
{
   // grav_vector.setVal(0.0);
    quantummoveKick(timestep, a_new, a_half,false);
    ParticleContainer<N>::moveKick(grav_vector, level, timestep, a_new, a_half, start_comp_for_accel);
    particles_moved = true;

}

template <int N>
void AxionParticleContainer<N>::quantummoveKick(Real timestep,Real a_newer,Real a_older,bool before_drift)
{
    const Geometry& geom = this->m_amr->Geom(0);
    const Real* dx       = geom.CellSize();
    const Real* plo      = geom.ProbLo();
    const Real half_dt   = 0.5 * timestep;
    const Real a_old_inv = 1 / a_older;
    Real max_acc         = 0.0;
    Real a_inv_sqrt = 1 / a_older / a_older;
    if(!before_drift)
      Real a_inv_sqrt = 1/ a_newer / a_newer;


    
    if (this->m_particles.size() == 0) return;
    PMap&      pmap      = this->m_particles[0]; 
    
    MultiFab pressure(this->m_amr->ParticleBoxArray(0), 1, smoothing_length);
    MultiFab pressuregrad(this->m_amr->ParticleBoxArray(0), 3, 2);
    
    computeQuantumPressure(pressuregrad,pressure,0);
    if(!geom.isAllPeriodic())
      pressure.setBndry(0.0);
    pressure.FillBoundary();
    geom.FillPeriodicBoundary(pressure);

    
    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = pressure[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id > 0)
            {
                BL_ASSERT(p.m_grid == grid);
                //
                // Note: m_data[0] is mass, 1 is v_x, ...
                //
                Real acc[BL_SPACEDIM];
        for(int is=0;is<BL_SPACEDIM;is++)
          acc[is] = 0.0;
                const IntVect lowcell(floor((p.m_pos[0]-plo[0])/dx[0] - smoothing_length) , floor((p.m_pos[1]-plo[1])/dx[1] - smoothing_length), floor((p.m_pos[2]-plo[2])/dx[2] - smoothing_length));
            Real roverh,x[3],h;
            for(int ks=0;ks<smoothing_length*2+1;ks++){
              for(int js=0;js<smoothing_length*2+1;js++){
           for(int is=0;is<smoothing_length*2+1;is++){
              IntVect cell = lowcell;
              cell[0] += is;
              cell[1] += js;
              cell[2] += ks;
              for(int d=0;d<BL_SPACEDIM;d++)
                        x[d] = -p.m_pos[d] + plo[d]+(cell[d]+ 0.5)*dx[d];
              h = (smoothing_length+0.5)*dx[0];
              roverh = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/h;
              if(roverh<=0.5){
                for(int d = 0;d<BL_SPACEDIM;d++){
              acc[d] += gfab(cell,0) * 8.0/M_PI/pow(h,3)*(-12.0*x[d]/h/h + 18.0*x[d]*roverh/h/h);
            //  std::cout << acc[d] << std::endl;
            }
              } else if (roverh<1.0) {
            for(int d = 0;d<BL_SPACEDIM;d++){
              acc[d] += gfab(cell,0) * 8.0/M_PI/pow(h,3)*-6.0*x[d]*(1.0-roverh)*(1.0-roverh)/roverh/h/h;
            }
              } 
           }
          }
        }
        
        for(int d=0;d<BL_SPACEDIM;d++){
          acc[d]*= dx[0]*dx[1]*dx[2];
          if(acc[d] > max_acc) max_acc = acc[d];
        }
                //
                // Define (a u)^new = (a u)^half + dt/2 grav^new
                //
                D_TERM(p.m_data[1] *= a_older;,
                       p.m_data[2] *= a_older;,
                       p.m_data[3] *= a_older;);

                D_TERM(p.m_data[1] += a_inv_sqrt * half_dt * acc[0] - damping_constant*half_dt*p.m_data[1];,
                       p.m_data[2] += a_inv_sqrt * half_dt * acc[1] - damping_constant*half_dt*p.m_data[2];,
                       p.m_data[3] += a_inv_sqrt * half_dt * acc[2] - damping_constant*half_dt*p.m_data[3];);

                //We have to divide by the older a here because the gravitational kick is yet to come
                D_TERM(p.m_data[1] *= a_old_inv;,
                       p.m_data[2] *= a_old_inv;,
                       p.m_data[3] *= a_old_inv;);
            }
        }
    }
    ParallelDescriptor::ReduceRealMax(max_acc);
    max_acceleration = max_acc;
    
}

template <int N>
Real AxionParticleContainer<N>::estTimestep (MultiFab& grav_vector, Real a, int lev, Real cfl)
{
    const Geometry& geom = this->m_amr->Geom(0);
    const Real* dx       = geom.CellSize();
    Real dt = 1e100;//1/sqrt(max_acceleration/dx[0]);
    Real dt2 = ParticleContainer<N>::estTimestep(grav_vector, a, lev, cfl); 
    if(dt <= dt2){
      if(ParallelDescriptor::IOProcessor()){
    std::cout << "timestep was determined by acceleration criterion" << std::endl;
      }
      return dt;
    }else{
      if(ParallelDescriptor::IOProcessor()){
    std::cout << "timestep was determined by velocity/gravity criterion" << std::endl;
      }
      return dt2;
    }
}

template <int N>
void AxionParticleContainer<N>::AssignDensity (PArray<MultiFab>& mf, int lev_min, int ncomp, int finest_level, double masstreshold) const
{
    ParticleContainer<N>::AssignDensity(mf,lev_min,ncomp,finest_level,masstreshold);
}

template <int N>
void
AxionParticleContainer<N>::AssignDensitySingleLevel (MultiFab& mf_to_be_filled,
                                                int       lev,
                                                int       ncomp,
                                                int       particle_lvl_offset) const
{

    BL_ASSERT(N >= 1);
    BL_ASSERT(ncomp == 1 || ncomp == BL_SPACEDIM+1);
    
    const Real      strttime    = ParallelDescriptor::second();
    
    if(particles_moved || ncomp > 1){      
    const Amr*      m_amr               = this->m_amr;
    const Array<PMap> m_particles       = this->m_particles;
    
    //MultiFab* mf_pointer;
    if(!initzd){
     this->mf_pointer = new MultiFab(m_amr->ParticleBoxArray(lev),BL_SPACEDIM+1,smoothing_length,Fab_allocate);
     initzd = true;
    }

    if (lev >= m_particles.size())
    {
        //
        // Don't do anything if there are no particles at this level.
        //
    if (mf_pointer != &mf_to_be_filled) delete mf_pointer;
        return;
    }

   // const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_amr->Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx_particle = m_amr->Geom(lev + particle_lvl_offset).CellSize();
    const Real*     dx          = gm.CellSize();
    const PMap&     pmap        = m_particles[lev];
    const int       n           = pmap.size();

    const int       M           = (2*smoothing_length+1)*(2*smoothing_length+1)*(2*smoothing_length+1);

    if (gm.isAnyPeriodic() && !gm.isAllPeriodic())
        amrex::Error("AssignDensity: problem must be periodic in no or all directions");

    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)
        (*mf_pointer)[mfi].setVal(0);
    //
    // This is a little funky.  What in effect this'll do is force
    // each thread to work on a single (separate) grid at a time.  That
    // way no thread will step on any other.  If there's only one grid per CPU,
    // then oh well ....
    //

    Array<int>         pgrd(n);
    Array<const PBox*> pbxs(n);

    int j = 0;
    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
         pmap_it != pmapEnd;
         ++pmap_it, ++j)
    {
        pgrd[j] =   pmap_it->first;
        pbxs[j] = &(pmap_it->second);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) if (n > 1)
#endif
    for (int i = 0; i < n; i++)
    {
        const PBox& pbx = *pbxs[i];
        FArrayBox&  fab = (*mf_pointer)[pgrd[i]];

        Array<Real>    fracs;
        Array<IntVect> cells;
        for (typename PBox::const_iterator it = pbx.begin(), End = pbx.end();
             it != End;
             ++it)
        {
            const ParticleType& p = *it;

            if (p.m_id <= 0) continue;

            fracs.resize(M);
            cells.resize(M);
         //   const IntVect lowcell(floor((p.m_pos[0]-plo[0])/dx[0] - smoothing_length + 0.5) , floor((p.m_pos[1]-plo[1])/dx[1] - smoothing_length + 0.5), floor((p.m_pos[2]-plo[2])/dx[2] - smoothing_length + 0.5));
        const IntVect lowcell(floor((p.m_pos[0]-plo[0])/dx[0] - smoothing_length) , floor((p.m_pos[1]-plo[1])/dx[1] - smoothing_length), floor((p.m_pos[2]-plo[2])/dx[2]- smoothing_length ));
        Real roverh;
//      for(int ks=0;ks<1;ks++){//smoothing_length*2;ks++){
//        for(int js=0;js<1;js++){//smoothing_length*2;js++){
//      for(int is=0;is<1;is++){//smoothing_length*2;is++){
        for(int ks=0;ks<smoothing_length*2+1;ks++){
          for(int js=0;js<smoothing_length*2+1;js++){
        for(int is=0;is<smoothing_length*2+1;is++){
          IntVect cell = lowcell;
          cell[0] += is;
          cell[1] += js;
          cell[2] += ks;
          
          cells[(smoothing_length*2+1)*(smoothing_length*2+1)*ks + (smoothing_length*2+1)*js + is] = cell;
          
          roverh = sqrt( (p.m_pos[1] - plo[1]-(cell[1]+ 0.5)*dx[1])*(p.m_pos[1] - plo[1]-(cell[1]+0.5)*dx[1]) + (p.m_pos[0] - plo[0]-(cell[0]+ 0.5)*dx[0])*(p.m_pos[0] - plo[0]-(cell[0]+0.5)*dx[0]) + (p.m_pos[2] - plo[2]-(cell[2]+ 0.5)*dx[2])*(p.m_pos[2] - plo[2]-(cell[2]+0.5)*dx[2]))/(smoothing_length+0.5)/dx[0];
         // roverh = p.m_pos[1];// - plo[0]-(cell[0]+ 0.5)*dx[0];
          if(roverh <= 0.5){
            fracs[(smoothing_length*2+1)*(smoothing_length*2+1)*ks + (smoothing_length*2+1)*js + is] = 8.0/M_PI/pow((smoothing_length+0.5)*dx[0],3)*(1.0-6.0*roverh*roverh+ 6.0 *roverh*roverh*roverh);
          }else if(roverh < 1.0){
            fracs[(smoothing_length*2+1)*(smoothing_length*2+1)*ks + (smoothing_length*2+1)*js + is] = 8.0/M_PI/pow((smoothing_length+0.5)*dx[0],3)*2.0*(1.0-roverh)*(1.0-roverh)*(1.0-roverh);
          }
        }
          }
        }
            //
            // If this is not fully periodic then we have to be careful that no
            // particle's support leaves the domain. We test this by checking the low
            // and high corners respectively.
            //
            //if (!gm.isAllPeriodic() && !this->allow_particles_near_boundary)
            //    if (!gm.Domain().contains(cells[0]) || !gm.Domain().contains(cells[M-1]))
            //        BoxLib::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");

            for (int i = 0; i < M; i++)
            {
                if (!fab.box().contains(cells[i])) continue;

                // If the domain is not periodic and we want to let particles
                //    live near the boundary but "throw away" the contribution that 
                //    does not fall into the domain ...
                if (!gm.isAllPeriodic() && this->allow_particles_near_boundary 
                                        && !gm.Domain().contains(cells[i])) continue;
                //
                // Sum up mass in first component.
                //
                fab(cells[i],0) += p.m_data[0] * fracs[i];
                //
                // Sum up momenta in next components.
                //
                for (int n = 1; n < ncomp; n++)
                   fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
            }
        }
    }

    mf_pointer->SumBoundary();
    gm.SumPeriodicBoundary(*mf_pointer);
    //
    // If ncomp > 1, first divide the momenta (component n) 
    // by the mass (component 0) in order to get velocities.
    // Be careful not to divide by zero.
    //
    for (int n = 1; n < ncomp; n++)
    {
        for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)
        {
            (*mf_pointer)[mfi].protected_divide((*mf_pointer)[mfi],0,n,1);
        }
    }
    particles_moved = false;
    }
    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
//  delete mf_pointer;
    }

    if (this->m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AssignDensity(single-level) time: " << stoptime << '\n';
        }
    }
}

template <int N>
void
AxionParticleContainer<N>::AssignDensityGradientsSingleLevel (MultiFab& mf_to_be_filled,
                                                int       lev) const
{

    BL_ASSERT(N >= 1);
    BL_ASSERT(mf_to_be_filled.nComp() == 3);
    BL_ASSERT(mf_to_be_filled.nGrow() == 4);
    
    const Amr*      m_amr               = this->m_amr;
    const Array<PMap> m_particles       = this->m_particles;
    
    MultiFab* mf_pointer;

    mf_pointer = &mf_to_be_filled;

    if (lev >= m_particles.size())
    {
        //
        // Don't do anything if there are no particles at this level.
        //
    if (mf_pointer != &mf_to_be_filled) delete mf_pointer;
        return;
    }

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_amr->Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx_particle = m_amr->Geom(lev).CellSize();
    const Real*     dx          = gm.CellSize();
    const PMap&     pmap        = m_particles[lev];
    const int       n           = pmap.size();
    
    const int       M           = 8*8*8;

    if (gm.isAnyPeriodic() && !gm.isAllPeriodic())
        amrex::Error("AssignDensity: problem must be periodic in no or all directions");

    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)
        (*mf_pointer)[mfi].setVal(0);
    //
    // This is a little funky.  What in effect this'll do is force
    // each thread to work on a single (separate) grid at a time.  That
    // way no thread will step on any other.  If there's only one grid per CPU,
    // then oh well ....
    //

    Array<int>         pgrd(n);
    Array<const PBox*> pbxs(n);

    int j = 0;
    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
         pmap_it != pmapEnd;
         ++pmap_it, ++j)
    {
        pgrd[j] =   pmap_it->first;
        pbxs[j] = &(pmap_it->second);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) if (n > 1)
#endif
    for (int i = 0; i < n; i++)
    {
        const PBox& pbx = *pbxs[i];
        FArrayBox&  fab = (*mf_pointer)[pgrd[i]];

        Array<Real>    fracsx,fracsy,fracsz;
        Array<IntVect> cells;
        for (typename PBox::const_iterator it = pbx.begin(), End = pbx.end();
             it != End;
             ++it)
        {
            const ParticleType& p = *it;

            if (p.m_id <= 0) continue;

            fracsx.resize(M);
            fracsy.resize(M);
            fracsz.resize(M);
            cells.resize(M);
            const IntVect lowcell(floor((p.m_pos[0]-plo[0])/dx[0] - 3.5) , floor((p.m_pos[1]-plo[1])/dx[1] - 3.5), floor((p.m_pos[2]-plo[2])/dx[2] - 3.5));
        Real roverh,x,y,z,h;
        for(int ks=0;ks<8;ks++){
          for(int js=0;js<8;js++){
        for(int is=0;is<8;is++){
          IntVect cell = lowcell;
          cell[0] += is;
          cell[1] += js;
          cell[2] += ks;
          
          cells[64*ks + 8*js + is] = cell;
          
          x = -p.m_pos[0] + plo[0]+(cell[0]+ 0.5)*dx[0];
          y = -p.m_pos[1] + plo[1]+(cell[1]+ 0.5)*dx[1];
          z = -p.m_pos[2] + plo[2]+(cell[2]+ 0.5)*dx[2];
          h = 4.0*dx[0];
          
          roverh = sqrt(x*x+y*y+z*z)/h;
          
          if(roverh <= 0.5){
            fracsx[64*ks + 8*js + is] = 8.0/M_PI/pow(h,3)*(-12.0*x/h/h + 18.0*x*roverh/h/h);
            fracsy[64*ks + 8*js + is] = 8.0/M_PI/pow(h,3)*(-12.0*y/h/h + 18.0*y*roverh/h/h);
            fracsz[64*ks + 8*js + is] = 8.0/M_PI/pow(h,3)*(-12.0*z/h/h + 18.0*z*roverh/h/h);
          }else if(roverh <= 1.0){
            fracsx[64*ks + 8*js + is] = 8.0/M_PI/pow(h,3)*-6.0*x*(1.0-roverh)*(1.0-roverh)/roverh/h/h;
            fracsy[64*ks + 8*js + is] = 8.0/M_PI/pow(h,3)*-6.0*y*(1.0-roverh)*(1.0-roverh)/roverh/h/h;
                fracsz[64*ks + 8*js + is] = 8.0/M_PI/pow(h,3)*-6.0*z*(1.0-roverh)*(1.0-roverh)/roverh/h/h;
          }
        }
          }
        }
            //
            // If this is not fully periodic then we have to be careful that no
            // particle's support leaves the domain. We test this by checking the low
            // and high corners respectively.
            //
            if (!gm.isAllPeriodic() && !this->allow_particles_near_boundary)
                if (!gm.Domain().contains(cells[0]) || !gm.Domain().contains(cells[M-1]))
                    amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");

            for (int i = 0; i < M; i++)
            {
                if (!fab.box().contains(cells[i])) continue;

                // If the domain is not periodic and we want to let particles
                //    live near the boundary but "throw away" the contribution that 
                //    does not fall into the domain ...
                if (!gm.isAllPeriodic() && this->allow_particles_near_boundary 
                                        && !gm.Domain().contains(cells[i])) continue;

                fab(cells[i],0) += p.m_data[0] * fracsx[i];
        fab(cells[i],1) += p.m_data[0] * fracsy[i];
                fab(cells[i],2) += p.m_data[0] * fracsz[i];
            }
        }
    }

    mf_pointer->SumBoundary();
    gm.SumPeriodicBoundary(*mf_pointer);

    if (this->m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AssignDensityGradients(single-level) time: " << stoptime << '\n';
        }
    }
}

template <int N>
void AxionParticleContainer<N>::doFiltering (MultiFab& mf,int comp, int niter, int nstride) const
{
  
    const Geometry& geom = this->m_amr->Geom(0);
    const Real* dx       = geom.CellSize();
       
    MultiFab largeDens(this->m_amr->ParticleBoxArray(0), 1, 4);
    largeDens.copy(mf,comp,0,1);
    largeDens.FillBoundary();
    geom.FillPeriodicBoundary(largeDens);
    
    double weights[4] = {1.0/8.0, 1.0/16.0, 1.0/32.0, 1.0/64.0};

    for (int stride=1;stride<=nstride;stride++){
        for (int iter=0;iter<niter;iter++){  
            for (MFIter mfi(largeDens); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.validbox();
                const int*  lo = box.loVect();
                const int*  hi = box.hiVect();
                int  ngrow = largeDens.nGrow();


                BL_FORT_PROC_CALL(LOWPASSFILTER, lowpassfilter)
                                  (BL_TO_FORTRAN(largeDens[mfi]),
                                  dx, lo, hi,&stride,&ngrow, weights);
            }
           
            largeDens.FillBoundary();
            geom.FillPeriodicBoundary(largeDens);
       }
    }
    
    mf.copy(largeDens,0,comp,1);
}

template <int N>
void AxionParticleContainer<N>::addQuantumPressure (MultiFab & acc_vector, int lev, Real a_old, Real a_half)
{
    // we need one ghostzone for cic...
    MultiFab axPressForceMf(this->m_amr->ParticleBoxArray(lev), 3, 2);
    MultiFab pressMf(this->m_amr->ParticleBoxArray(lev), 1, 3);
    this->computeQuantumPressure(axPressForceMf, pressMf, lev);

#ifdef AXPARTNOGRAV
    acc_vector.setVal(0.);
#endif
    
    //add acceleration to acc_vector
    MultiFab::Add(acc_vector, axPressForceMf, 0, 0, 3, 1);
}

template <int N>
void AxionParticleContainer<N>::addArtificialViscosity (MultiFab & acc_vector, int lev, Real a_old, Real a_half)
{
    const Geometry& geom = this->m_amr->Geom(0);
    const Real* dx       = geom.CellSize();
    const Real* plo      = geom.ProbLo();
    
    MultiFab pressureMf(this->m_amr->ParticleBoxArray(0), 3, 1);
    
    //Assign density and velocities
    MultiFab mf(this->m_amr->ParticleBoxArray(0), 4, 3);
    ParticleContainer<N>::AssignDensitySingleLevel(mf, 0, 4);
    mf.FillBoundary();
    geom.FillPeriodicBoundary(mf);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int*  lo = box.loVect();
        const int*  hi = box.hiVect();

        //calculate the pressure 
    //calculate the acceleration from the pressure
        BL_FORT_PROC_CALL(VISCOSITYPRESSURE, viscositypressure)
                      (BL_TO_FORTRAN(mf[mfi]), BL_TO_FORTRAN(pressureMf[mfi]),
               dx, lo, hi);
    }    
    //add acceleration to acc_vector
    MultiFab::Add(acc_vector, pressureMf, 0, 0, 3, 1);
}

template <int N>
void AxionParticleContainer<N>::computeQuantumPressure (MultiFab & axPressForceMf, MultiFab & pressureMf, int lev)
{
    const Geometry& geom = this->m_amr->Geom(0);
    const Real* dx       = geom.CellSize();
    
    //calculate density (call some funciton from original particle container)
    PArray<MultiFab> PMF;
    this->AssignDensity(PMF);
    //FIXME possibly average down (we only do single level for the moment)

    //largeDens contains the filtered density
    MultiFab largeDens(this->m_amr->ParticleBoxArray(lev), 1, 3);
    largeDens.copy(PMF[lev],0,0,1);
    if(!geom.isAllPeriodic())
      largeDens.setBndry(0.0);
    largeDens.FillBoundary();
    geom.FillPeriodicBoundary(largeDens);


    for (MFIter mfi(largeDens); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int*  lo = box.loVect();
        const int*  hi = box.hiVect();

        //calculate the pressure 
    //calculate the acceleration from the pressure
        BL_FORT_PROC_CALL(AXIONPRESSURE, axionpressure)
                      (BL_TO_FORTRAN(largeDens[mfi]), BL_TO_FORTRAN(axPressForceMf[mfi]),BL_TO_FORTRAN(pressureMf[mfi]),
               dx, lo, hi);
    }
    
#ifndef NDEBUG    
    if(axPressForceMf.contains_nan()){
        if(ParallelDescriptor::IOProcessor()){
        std::cout << "Pressure gradient contains Nan!" << std::endl;
    }
    }
    if(pressureMf.contains_nan()){
        if(ParallelDescriptor::IOProcessor()){
        std::cout << "Pressure contains Nan!" << std::endl;
    }
    }
#endif

}
    
template <int N>
void
AxionParticleContainer<N>::InitVarMass (MultiFab& mf)
{
    const int       MyProc      = ParallelDescriptor::MyProc();
    const int       IOProc      = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& geom        = this->m_amr->Geom(0);
    const Real*     dx          = geom.CellSize();

    
    this->m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    this->m_particles.resize(this->m_amr->finestLevel()+1);

    for (int lev = 0; lev < this->m_particles.size(); lev++)
    {
        BL_ASSERT(this->m_particles[lev].empty());
    }


    //
    // The grid should be initialized according to the ics...
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();
        const int   fab_ix = fab_hi[0] - fab_lo[0] + 1;
        const int   fab_jx = fab_hi[1] - fab_lo[1] + 1;
        const int   fab_kx = fab_hi[2] - fab_lo[2] + 1;

        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
            for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
                    ParticleType p;

                    IntVect indices(D_DECL(ix, jx, kx));

                for (int n = 0; n < BL_SPACEDIM; n++)
                {
                        //
            // A homogeneous distribution (for 1 p per cell in the center of the cell),
                        //
                    p.m_pos[n] = geom.ProbLo(n) + 
                                      (indices[n]+0.5)*dx[n];
                }
                
         //       p.m_pos[2] += 0.1 * dx[2];
        //    p.m_pos[0] += 0.1 * dx[2];

            // set mass from density
                p.m_data[0] = myFab(indices,0) * dx[0] * dx[1] * dx[2];
            if(p.m_data[0] < 1.0e-10) continue;

            // set velocity
                for (int n = 0; n < BL_SPACEDIM; n++)
                {
            p.m_data[n+1] = myFab(indices, n+1);
            }

                p.m_id      = ParticleBase::NextID();
                p.m_cpu     = MyProc;
    
                if (!ParticleBase::Where(p,this->m_amr))
                    {
                    ParticleBase::PeriodicShift(p,this->m_amr);

                        if (!ParticleBase::Where(p,this->m_amr))
                            amrex::Abort("ParticleContainer<N>::InitVarMass(): invalid particle");
            }

                BL_ASSERT(p.m_lev >= 0 && p.m_lev <= this->m_amr->finestLevel());
                //
                // Add it to the appropriate PBox at the appropriate level.
                //
                this->m_particles[p.m_lev][p.m_grid].push_back(p);
                }
            }
        }
    }
}

template <int N>
void
AxionParticleContainer<N>::InitVarCount (MultiFab& mf, long n_axpart)
{
    const int       MyProc      = ParallelDescriptor::MyProc();
    const int       IOProc      = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime    = ParallelDescriptor::second();
    const Amr*      M_amr       = this->m_amr;
    const Geometry& geom        = M_amr->Geom(0);
    const Real*     dx          = geom.CellSize();
    Array<PMap>&    M_particles = this->m_particles;
    int             npart;
    Real            r;
    
    Real            factor      = mf.norm1(0,0)/n_axpart; //compute density per particle

    this->m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    this->m_particles.resize(M_amr->finestLevel()+1);

    for (int lev = 0; lev < M_particles.size(); lev++)
    {
        BL_ASSERT(M_particles[lev].empty());
    }

    amrex::InitRandom(MyProc);
    //
    // The grid should be initialized according to the ics...
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();
        const int   fab_ix = fab_hi[0] - fab_lo[0] + 1;
        const int   fab_jx = fab_hi[1] - fab_lo[1] + 1;
        const int   fab_kx = fab_hi[2] - fab_lo[2] + 1;

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
                            ParticleType p;



                            for (int n = 0; n < BL_SPACEDIM; n++)
                            {
                                do
                                {
                                    r = amrex::Random();
                                }
                                while (r == 0 || r == 1);

                                p.m_pos[n] = geom.ProbLo(n) +
                                    (indices[n]+r)*dx[n];
                            }

                            // set mass
                            p.m_data[0] = factor * dx[0] * dx[1] * dx[2];

                            // set velocity
                            for (int n = 0; n < BL_SPACEDIM; n++)
                            {
                                p.m_data[n+1] = myFab(indices, n+1);
                            }

                            p.m_id      = ParticleBase::NextID();
                            p.m_cpu     = MyProc;

                            if (!ParticleBase::Where(p,M_amr))
                            {
                                ParticleBase::PeriodicShift(p,M_amr);

                            if (!ParticleBase::Where(p,M_amr))
                                amrex::Abort("ParticleContainer<N>::InitVarCount(): invalid particle");
                            }

                            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= M_amr->finestLevel());
                            //
                            // Add it to the appropriate PBox at the appropriate level.
                            //
                            M_particles[p.m_lev][p.m_grid].push_back(p);
             }
                }
            }
        }
    }
}

#endif //defined _AXPARTICLES_H_


````

