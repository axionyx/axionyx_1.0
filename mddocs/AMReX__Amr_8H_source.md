
# File AMReX\_Amr.H

[**File List**](files.md) **>** [**AMReX\_axionyx**](dir_5c77c3c750fcf9b051dca9dbb6924de0.md) **>** [**AMReX\_Amr.H**](AMReX__Amr_8H.md)

[Go to the documentation of this file.](AMReX__Amr_8H.md) 


````cpp

#ifndef AMREX_Amr_H_
#define AMREX_Amr_H_

#include <fstream>
#include <memory>
#include <list>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BCRec.H>

#include <AMReX_AmrCore.H>

#ifdef USE_PERILLA
#include <RegionGraph.H>
#include <Perilla.H>
#endif

namespace amrex {

class AmrLevel;
class LevelBld;
class BoxDomain;
template <class T>
class MFGraph;
class AmrTask;
#if defined(BL_USE_SENSEI_INSITU)
class AmrInSituBridge;
#endif

class Amr
    : public AmrCore
{
  template <class T>
  friend class MFGraph;
  friend class AmrTask;
  typedef std::multimap< std::pair<int, int>, double >  BoundaryPointList;

public:
    Amr ();

    Amr (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord);

    Amr (const Amr& rhs) = delete;
    Amr& operator= (const Amr& rhs) = delete;

    void InitAmr ();

    virtual ~Amr ();

    virtual void init (Real strt_time, Real stop_time);

    void InitializeInit (Real strt_time, Real stop_time,
                         const BoxArray* lev0_grids = 0, const Vector<int>* pmap = 0);

    void FinalizeInit (Real strt_time, Real stop_time);

    void setDtLevel (const Vector<Real>& dt_lev) noexcept;

    void setDtLevel (Real dt, int lev) noexcept;

    void setDtMin (const Vector<Real>& dt_lev) noexcept;

    void setNCycle (const Vector<int>& mss) noexcept;

    int subCycle () const noexcept { return sub_cycle; }

    const std::string& subcyclingMode() const noexcept { return subcycling_mode; }

    int level_being_advanced () const noexcept { return which_level_being_advanced; }
    Real cumTime () const noexcept { return cumtime; }
    void setCumTime (Real t) noexcept {cumtime = t;}
    Real startTime () const noexcept { return start_time; }
    void setStartTime (Real t) noexcept {start_time = t;}
    Real dtLevel (int level) const noexcept { return dt_level[level]; }
    Real dtMin (int level) const noexcept { return dt_min[level]; }
    const Vector<Real>& dtLevel () const noexcept { return dt_level; }
    int nCycle (int level) const noexcept { return n_cycle[level]; }
    int levelSteps (int lev) const noexcept { return level_steps[lev]; }
    void setLevelSteps (int lev, int n) noexcept { level_steps[lev] = n; }
    int levelCount (int lev) const noexcept { return level_count[lev]; }
    void setLevelCount (int lev, int n) noexcept { level_count[lev] = n; }
    bool RegridOnRestart () const noexcept;
    int regridInt (int lev) const noexcept { return regrid_int[lev]; }
    int checkInt () const noexcept { return check_int; }
    Real checkPer() const noexcept { return check_per; }
    int plotInt () const noexcept { return plot_int; }
    Real plotPer () const noexcept { return plot_per; }
    Real plotLogPer () const noexcept { return plot_log_per; }
    int smallplotInt () const noexcept { return small_plot_int; }
    Real smallplotPer () const noexcept { return small_plot_per; }
    Real smallplotLogPer () const noexcept { return small_plot_log_per; }
    static const std::list<std::string>& statePlotVars () noexcept { return state_plot_vars; }
    static const std::list<std::string>& stateSmallPlotVars () noexcept { return state_small_plot_vars; }
    static bool isStatePlotVar (const std::string& name);
    static bool isStateSmallPlotVar (const std::string& name);
    static void addStatePlotVar (const std::string& name);
    static void addStateSmallPlotVar (const std::string& name);
    static void deleteStatePlotVar (const std::string& name);
    static void clearStatePlotVarList ();
    static void clearStateSmallPlotVarList ();
    static void fillStatePlotVarList ();
    static void fillStateSmallPlotVarList ();
    static bool Plot_Files_Output ();
    static const std::list<std::string>& derivePlotVars () noexcept { return derive_plot_vars; }
    static const std::list<std::string>& deriveSmallPlotVars () noexcept { return derive_small_plot_vars; }
    static bool isDerivePlotVar (const std::string& name) noexcept;
    static bool isDeriveSmallPlotVar (const std::string& name) noexcept;
    static void addDerivePlotVar (const std::string& name);
    static void addDeriveSmallPlotVar (const std::string& name);
    static void deleteDerivePlotVar (const std::string& name);
    static void deleteDeriveSmallPlotVar (const std::string& name);
    static void clearDerivePlotVarList ();
    static void clearDeriveSmallPlotVarList ();
    static void fillDerivePlotVarList ();
    static void fillDeriveSmallPlotVarList ();

    static void Initialize ();
    static void Finalize ();
    AmrLevel& getLevel (int lev) noexcept { return *amr_level[lev]; }
    Vector<std::unique_ptr<AmrLevel> >& getAmrLevels () noexcept;
    long cellCount () noexcept;
    long cellCount (int lev) noexcept;
    int numGrids () noexcept;
    int numGrids (int lev) noexcept;
    int okToContinue () noexcept;
    void RegridOnly (Real time, bool do_io = true);
    bool okToRegrid (int level) noexcept;
    static const BoxArray& initialBa (int level) noexcept
        { BL_ASSERT(level-1 < initial_ba.size()); return initial_ba[level-1]; }
    static int initialBaLevels () noexcept { return initial_ba.size(); }
    virtual void coarseTimeStep (Real stop_time);

    Real coarseTimeStepDt (Real stop_time);
    std::unique_ptr<MultiFab> derive (const std::string& name,
                      Real           time,
                      int            lev,
                      int            ngrow);
    const std::string& theRestartFile () const noexcept { return restart_chkfile; }
    const std::string& theRestartPlotFile () const noexcept { return restart_pltfile; }
    std::ostream& DataLog (int i);
    const std::string DataLogName (int i) const noexcept { return datalogname[i]; }
    int NumDataLogs () noexcept;
    static Real computeOptimalSubcycling (int   n,
                                          int*  best,
                                          Real* dt_max,
                                          Real* est_work,
                                          int*  cycle_max);

    virtual void writePlotFile ();
    int stepOfLastPlotFile () const noexcept {return last_plotfile;}
    virtual void writeSmallPlotFile ();
    int stepOfLastSmallPlotFile () const noexcept {return last_smallplotfile;}
    virtual void checkPoint ();
    int stepOfLastCheckPoint () const noexcept {return last_checkpoint;}

    const Vector<BoxArray>& getInitialBA() noexcept;

    void setBoundaryGeometry(BoundaryPointList& IntersectLoX,
                             BoundaryPointList& IntersectHiX,
                             BoundaryPointList& IntersectLoY,
                             BoundaryPointList& IntersectHiY) noexcept
    {
        intersect_lox = IntersectLoX;
        intersect_hix = IntersectHiX;
        intersect_loy = IntersectLoY;
        intersect_hiy = IntersectHiY;
    };

    void setBoundaryGeometry(BoundaryPointList& IntersectLoX,
                             BoundaryPointList& IntersectHiX,
                             BoundaryPointList& IntersectLoY,
                             BoundaryPointList& IntersectHiY,
                             BoundaryPointList& IntersectLoZ,
                             BoundaryPointList& IntersectHiZ) noexcept
    {
        intersect_lox = IntersectLoX;
        intersect_hix = IntersectHiX;
        intersect_loy = IntersectLoY;
        intersect_hiy = IntersectHiY;
        intersect_loz = IntersectLoZ;
        intersect_hiz = IntersectHiZ;
    };

    BoundaryPointList& getIntersectLoX() noexcept
    {
        return intersect_lox;
    };
    BoundaryPointList& getIntersectHiX() noexcept
    {
        return intersect_hix;
    };
    BoundaryPointList& getIntersectLoY() noexcept
    {
        return intersect_loy;
    };
    BoundaryPointList& getIntersectHiY() noexcept
    {
        return intersect_hiy;
    };
    BoundaryPointList& getIntersectLoZ() noexcept
    {
        return intersect_loz;
    };
    BoundaryPointList& getIntersectHiZ() noexcept
    {
        return intersect_hiz;
    };

#ifdef AMREX_PARTICLES
    void RedistributeParticles ();
#endif

    void InstallNewDistributionMap (int lev, const DistributionMapping& newdm);

    bool UsingPrecreateDirectories () noexcept;

protected:

    void initialInit (Real strt_time, Real stop_time,
                      const BoxArray* lev0_grids = 0, const Vector<int>* pmap = 0);
    void readProbinFile (int& init);
    void checkInput ();
    void restart (const std::string& filename);
    void defBaseLevel (Real start_time, const BoxArray* lev0_grids = 0, const Vector<int>* pmap = 0);
    void bldFineLevels (Real start_time);
    virtual void regrid (int  lbase,
             int  iteration,
                         Real time,
                         bool initial = false) override;
    virtual void regrid_level_0_on_restart ();
    void grid_places (int              lbase,
                      Real             time,
                      int&             new_finest,
                      Vector<BoxArray>& new_grids);

    DistributionMapping makeLoadBalanceDistributionMap (int lev, Real time, const BoxArray& ba) const;
    void LoadBalanceLevel0 (Real time);

    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override;
    virtual BoxArray GetAreaNotToTag (int lev) override;
    virtual void ManualTagsPlacement (int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev) override;

    virtual void timeStep (int  level,
                           Real time,
                           int  iteration,
                           int  niter,
                           Real stop_time);

    // pure virtural function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override
    { amrex::Abort("How did we get her!"); }
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override
    { amrex::Abort("How did we get her!"); }
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override
    { amrex::Abort("How did we get her!"); }
    virtual void ClearLevel (int lev) override
    { amrex::Abort("How did we get her!"); }

    bool writePlotNow () noexcept;
    bool writeSmallPlotNow () noexcept;

    void printGridInfo (std::ostream& os,
                        int           min_lev,
                        int           max_lev);

    void setRecordGridInfo (const std::string&);

    void setRecordRunInfo (const std::string&);

    void setRecordRunInfoTerse (const std::string&);

    void setRecordDataInfo (int i, const std::string&);

    void initSubcycle();
    void initPltAndChk();

    int initInSitu();
    int updateInSitu();
    static int finalizeInSitu();

    //
    // The data ...
    //
    std::string      regrid_grids_file;   
    std::string      initial_grids_file;  
    Vector<std::unique_ptr<AmrLevel> > amr_level;    
    Real             cumtime;      
    Real             start_time;   
    Vector<Real>      dt_level;     
    Vector<int>       level_steps;  
    Vector<int>       level_count;
    Vector<int>       n_cycle;
    std::string      subcycling_mode; 
    Vector<Real>      dt_min;
    bool             isPeriodic[AMREX_SPACEDIM];  
    Vector<int>       regrid_int;      
    int              last_checkpoint; 
    int              check_int;       
    Real             check_per;       
    std::string      check_file_root; 
    int              last_plotfile;   
    int              last_smallplotfile;   
    int              plot_int;        
    Real             plot_per;        
    Real             plot_log_per;    
    int              small_plot_int;  
    Real             small_plot_per;  
    Real             small_plot_log_per;  
    int              write_plotfile_with_checkpoint;  
    int              file_name_digits; 
    int              message_int;     
    std::string      plot_file_root;  
    std::string      small_plot_file_root;  

    int              which_level_being_advanced; 

    int              record_grid_info;
    int              record_run_info;
    int              record_run_info_terse;
    std::ofstream    gridlog;
    std::ofstream    runlog;
    std::ofstream    runlog_terse;
    Vector<std::unique_ptr<std::fstream> > datalog;
    Vector<std::string> datalogname;
    int              sub_cycle;
    std::string      restart_chkfile;
    std::string      restart_pltfile;
    std::string      probin_file;
    LevelBld*        levelbld;
    bool             abort_on_stream_retry_failure;
    int              stream_max_tries;
    int              loadbalance_with_workestimates;
    int              loadbalance_level0_int;
    Real             loadbalance_max_fac;

    bool             bUserStopRequest;

    //
    // The static data ...
    //
    static std::list<std::string> state_plot_vars;  
    static std::list<std::string> state_small_plot_vars;  
    static std::list<std::string> derive_plot_vars; 
    static std::list<std::string> derive_small_plot_vars; 
    static bool                   first_plotfile;
    static Vector<BoxArray> initial_ba;
    static Vector<BoxArray> regrid_ba;

#if defined(BL_USE_SENSEI_INSITU)
    static AmrInSituBridge *insitu_bridge;
#endif

public:
#ifdef USE_PERILLA
    std::vector<std::vector<RegionGraph*> > graphArray;
    std::vector<RegionGraph*> amrGraphArray;
    std::vector<RegionGraph*> &get_graphArray(int level){return graphArray[level];}
#endif

    BoundaryPointList intersect_lox;
    BoundaryPointList intersect_loy;
    BoundaryPointList intersect_loz;
    BoundaryPointList intersect_hix;
    BoundaryPointList intersect_hiy;
    BoundaryPointList intersect_hiz;

    static bool first_smallplotfile;

};

}

#endif /*_Amr_H_*/
````

