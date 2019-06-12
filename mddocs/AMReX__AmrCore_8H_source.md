
# File AMReX\_AmrCore.H

[**File List**](files.md) **>** [**AMReX\_axionyx**](dir_5c77c3c750fcf9b051dca9dbb6924de0.md) **>** [**AMReX\_AmrCore.H**](AMReX__AmrCore_8H.md)

[Go to the documentation of this file.](AMReX__AmrCore_8H.md) 


````cpp
#ifndef BL_AMRCORE_H_
#define BL_AMRCORE_H_

#include <ostream>
#include <memory>

#include <AMReX_AmrMesh.H>

namespace amrex {

#ifdef AMREX_PARTICLES
class AmrParGDB;
#endif

class AmrCore
    : public AmrMesh
{
public:

    AmrCore ();
    AmrCore (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord=-1, Vector<IntVect> ref_ratios = Vector<IntVect>());

    AmrCore (const AmrCore& rhs) = delete;
    AmrCore& operator= (const AmrCore& rhs) = delete;

    virtual ~AmrCore ();

#ifdef AMREX_PARTICLES
    AmrParGDB* GetParGDB () const noexcept { return m_gdb.get(); }
#endif

    void InitFromScratch (Real time);

    virtual void regrid (int lbase, int iteration, Real time, bool initial=false);

    static void Initialize ();
    static void Finalize ();

    void printGridSummary (std::ostream& os, int min_lev, int max_lev) const noexcept;

    int Verbose () const noexcept { return verbose; }

protected:

    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override = 0;

    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override = 0;

    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    virtual void ClearLevel (int lev) = 0;

    int              verbose;

#ifdef AMREX_PARTICLES
    std::unique_ptr<AmrParGDB> m_gdb;
#endif

private:
    void InitAmrCore ();
};

}

#endif
````

