#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
using namespace amrex;
void divergence(Box const& bx, Array4<Real> const& div,
		Array4<Real const> const& vel, 
		const amrex::Real* dx, const amrex::Real a);
