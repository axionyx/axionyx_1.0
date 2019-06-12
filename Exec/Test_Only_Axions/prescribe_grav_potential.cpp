#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
using namespace amrex;
#ifdef CGRAV
void prescribe_grav_potential(MultiFab& phi,const Geometry &geom, int level, int finest_level)
{
    //once we have updated amrex, we can use the syntax below.
    const int our_comp = 0;
    const Real *ProbLo = geom.ProbLo();
    const Real *ProbHi = geom.ProbHi();
    const Real *dx = geom.CellSize();
    for (MFIter mfi(phi,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Dim3 lo = amrex::lbound(bx);
        const Dim3 hi = amrex::ubound(bx);
        Array4<Real> const& src = phi[mfi].array();
        for(int k = lo.z; k <= hi.z; ++k)
        {
            Real z = ProbLo[2]+(k+0.5)*dx[2];
            for(int j = lo.y; j <= hi.y; ++j)
            {
                Real y = ProbLo[1]+(j+0.5)*dx[1];
                AMREX_PRAGMA_SIMD
                for(int i = lo.x; i <= hi.x; ++i) 
                {
                    Real x = ProbLo[0]+(i+0.5)*dx[0];
                    //here, x, y, z are the physical coordinates of the cell (i,j,k).
                    Real &value = src(i,j,k,our_comp);
                    //TODO do something to value
                    
                    //for now, just set it to zero.
                    value = 0.0;

                }

            }

        }




    }

}
#endif
