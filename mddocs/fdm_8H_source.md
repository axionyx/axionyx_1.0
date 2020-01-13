
# File fdm.H

[**File List**](files.md) **>** [**FDM**](dir_43b815edcf2a06ee60d8a45cc6c77fb8.md) **>** [**fdm.H**](fdm_8H.md)

[Go to the documentation of this file.](fdm_8H.md) 


````cpp
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
using namespace amrex;
void divergence(Box const& bx, Array4<Real> const& div,
        Array4<Real const> const& vel, 
        const amrex::Real* dx, const amrex::Real a);
````

