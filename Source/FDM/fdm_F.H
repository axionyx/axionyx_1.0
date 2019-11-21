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

    void update_gaussian_beams(const int* np, void* particles,
			      const amrex::Real* accel, const int* accel_lo, const int* accel_hi,
			      const amrex::Real* phi, const int* phi_lo, const int* phi_hi,
			      const amrex::Real* prob_lo, 
			      const amrex::Real* dx, const amrex::Real& dt, 
			       const amrex::Real& a_prev, const amrex::Real& a_cur, const int* do_move);

    void update_fdm_particles_wkb(const int* np, void* particles,
				  const amrex::Real* accel, const int* accel_lo, const int* accel_hi,
				  const amrex::Real* prob_lo, 
				  const amrex::Real* dx, const amrex::Real& dt, 
				  const amrex::Real& a_prev, const amrex::Real& a_cur, const int* do_move);

    void update_gaussian_beams_wkb(const int* np, void* particles,
				   const amrex::Real* accel, const int* accel_lo, const int* accel_hi,
				   const amrex::Real* phi, const int* phi_lo, const int* phi_hi,
				   const amrex::Real* prob_lo, 
				   const amrex::Real* dx, const amrex::Real& dt, 
				   const amrex::Real& a_prev, const amrex::Real& a_cur, const int* do_move);

  void divergence( const amrex::Real* div, const amrex::Real* vel, const int* vel_lo, const int* vel_hi, const amrex::Real* dx);

  void initpart(const int* level, const double* time,
		const int* lo, const int* hi,
		const int* nd, BL_FORT_FAB_ARG(dat),
		const double* delta,
		const double* xlo, const double* xhi);

  void deposit_fdm_particles(const void* particles, const long* np, const amrex::Real* state_real,
			     const int* lo_real, const int* hi_real, const amrex::Real* state_imag,
			     const int* lo_imag, const int* hi_imag, const amrex::Real* prob_lo,
			     const amrex::Real* dx, const amrex::Real& a);

  void deposit_fdm_particles_wkb(const void* particles, const long* np, const amrex::Real* state_real,
				 const int* lo_real, const int* hi_real, const amrex::Real* state_imag,
				 const int* lo_imag, const int* hi_imag, const amrex::Real* prob_lo,
				 const amrex::Real* dx, const amrex::Real& a);

  void fort_set_mtt(const amrex::Real& m_tt);
  void fort_set_hbaroverm(const amrex::Real& hbaroverm);
  void fort_set_theta(const amrex::Real& theta_fdm);
  void fort_set_sigma(const amrex::Real& sigma_fdm);
  void fort_set_gamma(const amrex::Real& gamma_fdm);
  void fort_set_meandens(const amrex::Real& meandens);
  void fort_set_a(const amrex::Real& a);
  void fort_set_ratio(const amrex::Real& ratio_fdm);

#ifdef __cplusplus
}
#endif
