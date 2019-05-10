#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
using namespace amrex;
#ifdef CGRAV
void prescribe_grav_potential(MultiFab& phi, Geometry &geom, int level, int finest_level)
{
    //once we have updated amrex, we can use the syntax below.
//    const int our_comp = 0;
//    for (MFIter mfi(phi,true); mfi.isValid(); ++mfi)
//    {
//        const Box& bx = mfi.tilebox();
//        const Dim3 lo = amrex::lbound(bx);
//        const Dim3 hi = amrex::ubound(bx);
//        Array4<Real> const& src = phi[mfi].array();
//        for(int k = lo.z; k <= hi.z; ++k)
//        {
//            for(int j = lo.y; j <= hi.y; ++j)
//            {
//                for(int i = lo.x; i <= hi.x; ++i) 
//                {
//                    Real &value = src(i,j,k,our_comp);
//                    //TODO do something to value
//
//                }
//
//            }
//
//        }
//
//
//
//
//    }

}
#endif
