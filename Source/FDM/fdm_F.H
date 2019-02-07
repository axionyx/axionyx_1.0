#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void update_fdm_particles(const int* np, void* particles,
			      const amrex::Real* accel, const int* accel_lo, const int* accel_hi,
			      const amrex::Real* prob_lo, 
			      const amrex::Real* dx, const amrex::Real& dt, 
			      const amrex::Real& a_prev, const amrex::Real& a_cur, const int* do_move);

  void initpart(const int* level, const double* time,
		const int* lo, const int* hi,
		const int* nd, BL_FORT_FAB_ARG(dat),
		const double* delta,
		const double* xlo, const double* xhi);

  void deposit_fdm_particles(const void* particles, const long* np, const amrex::Real* state_real,
			     const int* lo_real, const int* hi_real, const amrex::Real* state_imag,
			     const int* lo_imag, const int* hi_imag, const double* xlo, const double* delta);
#ifdef __cplusplus
}
#endif
