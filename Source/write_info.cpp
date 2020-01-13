#include <iomanip>
#include "Nyx.H"

using namespace amrex;

void
Nyx::write_info ()
{
    int ndatalogs = parent->NumDataLogs();
    Real time_unit = 1.0;//3.0856776e19 / 31557600.0; // conversion to Julian years

    if (ndatalogs > 0)
    {
#ifndef NO_HYDRO
        MultiFab& S_new = get_new_data(State_Type);
        MultiFab& D_new = get_new_data(DiagEOS_Type);
	Real      max_t = 0;

        Real rho_T_avg=0.0, T_avg=0.0, Tinv_avg=0.0, T_meanrho=0.0;
        Real whim_mass_frac, whim_vol_frac, hh_mass_frac, hh_vol_frac, igm_mass_frac, igm_vol_frac;

	if (do_hydro)
        {
            // Removed reset internal energy before call to compute_temp, still compute new temp
    	    compute_new_temp(S_new,D_new);
            max_t = D_new.norm0(Temp_comp);
            compute_rho_temp(rho_T_avg, T_avg, Tinv_avg, T_meanrho);
            compute_gas_fractions(1.0e5, 120.0, whim_mass_frac, whim_vol_frac,
                                  hh_mass_frac, hh_vol_frac, igm_mass_frac, igm_vol_frac);
	}
#endif
#ifdef FDM
        Real mass=0.0, epot=0.0, ekinrho=0.0, ekinv=0.0, ekin=0.0, etot=0.0;
        Real angmom_x=0.0, angmom_y=0.0, angmom_z=0.0, grav_pot=0.0, phase=0.0;
        compute_axion_quantities(mass, epot, ekinrho, ekinv, ekin, angmom_x, angmom_y, angmom_z, grav_pot, phase);
	etot = epot + ekin;//ekinrho + ekinv;
	Real max_dens = 0.0, max_dens_nb = 0.0;
	for (int lev = 0; lev <= parent->finestLevel(); lev++)
	  max_dens = std::max(max_dens, get_level(lev).get_new_data(Axion_Type).max(Nyx::AxDens));
	for (int lev = 0; lev <= parent->finestLevel(); lev++)
	  max_dens_nb = std::max(max_dens_nb, get_level(lev).derive("particle_mass_density", state[Axion_Type].curTime(), 0)->max(0));
#endif

#ifdef NO_HYDRO
        Real time  = state[PhiGrav_Type].curTime();
#else
        Real time  = state[State_Type].curTime();
#endif
        Real dt    = parent->dtLevel(0);
        int  nstep = parent->levelSteps(0);

        if (ParallelDescriptor::IOProcessor())
        {
            std::ostream& data_loga = parent->DataLog(0);

            if (time == 0.0)
            {
                data_loga << std::setw( 8) <<  "#  nstep";
                data_loga << std::setw(14) <<  "       time   ";
                data_loga << std::setw(14) <<  "       dt     ";
                data_loga << std::setw(14) <<  "         z    ";
#ifdef FDM
                 data_loga << std::setw(14) <<  "     Mass";
                 data_loga << std::setw(14) <<  "     Epot";
                 data_loga << std::setw(14) <<  "  Ekinrho";
                 data_loga << std::setw(14) <<  "    Ekinv";
                 data_loga << std::setw(14) <<  "     Ekin";
                 data_loga << std::setw(14) <<  "     Etot";
                 data_loga << std::setw(22) <<  " max_dens";
                 data_loga << std::setw(22) <<  " max_dens_nb";
#endif
                 data_loga << std::setw(14) <<  "        a";
#ifndef NO_HYDRO
                if (do_hydro == 1)
                {
                   data_loga << std::setw(14) <<  "    T_max      ";
                   data_loga << std::setw(14) <<  "  <T>_rho      ";
                   data_loga << std::setw(14) <<  "  <T>_V        ";
                   data_loga << std::setw(14) <<  "T @ <rho>      ";
                   data_loga << std::setw(14) <<  "T(21cm)        ";
                   data_loga << std::setw(14) <<  "adiab.         ";
                   data_loga << std::setw(14) <<  "WHIM_m         ";
                   data_loga << std::setw(14) <<  "WHIM_v         ";
                   data_loga << std::setw(14) <<  "HH_m           ";
                   data_loga << std::setw(14) <<  "HH_v           ";
                   data_loga << std::setw(14) <<  "IGM_m          ";
                   data_loga << std::setw(14) <<  "IGM_v          ";
                }
#endif
                data_loga << '\n';

                Real old_z = (1. / old_a) - 1.;
                data_loga << std::setw( 8) <<  nstep;
                data_loga << std::setw(14) <<  std::setprecision(6) <<  (initial_time+time) * time_unit;
                data_loga << std::setw(14) <<  std::setprecision(6) <<    dt * time_unit;
                data_loga << std::setw(14) <<  std::setprecision(6) << old_z;
#ifdef FDM
                 data_loga << std::setw(14) <<  std::setprecision(6) << mass;
                 data_loga << std::setw(14) <<  std::setprecision(6) << epot;
                 data_loga << std::setw(14) <<  std::setprecision(6) << ekinrho;
                 data_loga << std::setw(14) <<  std::setprecision(6) << ekinv;
                 data_loga << std::setw(14) <<  std::setprecision(6) << ekin;
                 data_loga << std::setw(14) <<  std::setprecision(6) << etot;
                 data_loga << std::setw(22) <<  std::setprecision(6) << max_dens;
                 data_loga << std::setw(22) <<  std::setprecision(6) << max_dens_nb;
#endif
                data_loga << std::setw(14) <<  std::setprecision(6) << old_a;
#ifndef NO_HYDRO
                if (do_hydro == 1)
                {
                   data_loga << std::setw(14) <<  std::setprecision(6) << max_t;
                   data_loga << std::setw(14) <<  std::setprecision(6) << rho_T_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << T_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << T_meanrho;
                   data_loga << std::setw(14) <<  std::setprecision(6) << 1.0/Tinv_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << 0.021*(1.0+old_z)*(1.0+old_z);
                   data_loga << std::setw(14) <<  std::setprecision(6) << whim_mass_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << whim_vol_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << hh_mass_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << hh_vol_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << igm_mass_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << igm_vol_frac;
                }
#endif
                data_loga << '\n';
            }
            else
            {
                const Real new_z = (1. / new_a) - 1.;
                data_loga << std::setw( 8) <<  nstep;
                data_loga << std::setw(14) <<  std::setprecision(6) <<  (initial_time+time) * time_unit;
                data_loga << std::setw(14) <<  std::setprecision(6) <<    dt * time_unit;
                data_loga << std::setw(14) <<  std::setprecision(6) << new_z;

#ifdef FDM
                 data_loga << std::setw(14) <<  std::setprecision(6) << mass;
                 data_loga << std::setw(14) <<  std::setprecision(6) << epot;
                 data_loga << std::setw(14) <<  std::setprecision(6) << ekinrho;
                 data_loga << std::setw(14) <<  std::setprecision(6) << ekinv;
                 data_loga << std::setw(14) <<  std::setprecision(6) << ekin;
                 data_loga << std::setw(14) <<  std::setprecision(6) << etot;
                 data_loga << std::setw(22) <<  std::setprecision(10) << max_dens;
                 data_loga << std::setw(22) <<  std::setprecision(10) << max_dens_nb;
#endif
                data_loga << std::setw(14) <<  std::setprecision(6) << new_a;
#ifndef NO_HYDRO
                if (do_hydro == 1)
                {
                   data_loga << std::setw(14) <<  std::setprecision(6) << max_t;
                   data_loga << std::setw(14) <<  std::setprecision(6) << rho_T_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << T_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << T_meanrho;
                   data_loga << std::setw(14) <<  std::setprecision(6) << 1.0/Tinv_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << 0.021*(1.0+new_z)*(1.0+new_z);
                   data_loga << std::setw(14) <<  std::setprecision(6) << whim_mass_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << whim_vol_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << hh_mass_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << hh_vol_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << igm_mass_frac;
                   data_loga << std::setw(14) <<  std::setprecision(6) << igm_vol_frac;
                }
#endif
                data_loga << std::endl;
            }
        }
    }
}
