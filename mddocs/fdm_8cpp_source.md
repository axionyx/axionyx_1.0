
# File fdm.cpp

[**File List**](files.md) **>** [**FDM**](dir_43b815edcf2a06ee60d8a45cc6c77fb8.md) **>** [**fdm.cpp**](fdm_8cpp.md)

[Go to the documentation of this file.](fdm_8cpp.md) 


````cpp
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
using namespace amrex;
void divergence (Box const& bx, Array4<Real> const& div,
         Array4<Real const> const& vel, Real const* dx, const Real a)
{
  const auto lo = lbound(bx);
  const auto hi = ubound(bx);
  for     (int k = lo.z; k <= hi.z; ++k) {
    for   (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
    div(i,j,k,0) = (vel(i+1,j,k,0)-vel(i-1,j,k,0))/dx[0]*a/2.0
                  +(vel(i,j+1,k,1)-vel(i,j-1,k,1))/dx[1]*a/2.0
                  +(vel(i,j,k+1,2)-vel(i,j,k-1,2))/dx[2]*a/2.0;
      }
    }
  }
}
````

