module particle_mod

  use amrex_fort_module, only: c_real => amrex_particle_real
  use iso_c_binding ,    only: c_int

  implicit none
  private

  public  dm_particle_t
  
  type, bind(C)  :: dm_particle_t
     real(c_real)    :: pos(3)     !< Position
     real(c_real)    :: mass       !< Particle mass
     real(c_real)    :: vel(3)     !< Particle velocity
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type dm_particle_t

  public  agn_particle_t
  
  type, bind(C)  :: agn_particle_t
     real(c_real)    :: pos(3)     !< Position
     real(c_real)    :: mass       !< Particle mass
     real(c_real)    :: vel(3)     !< Particle velocity
     real(c_real)    :: energy     !< Particle energy
     real(c_real)    :: mdot       !< Particle mass change
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type agn_particle_t

  public  fdm_particle_t
  
  type, bind(C)  :: fdm_particle_t
     real(c_real)    :: pos(3)     !< Position
     real(c_real)    :: mass       !< Particle mass
     real(c_real)    :: vel(3)     !< Particle velocity
     real(c_real)    :: phase      !< Particle phase 
     real(c_real)    :: amp(2)     !< Particle amplitude
     real(c_real)    :: width      !< Particle width
     real(c_real)    :: qq(9)      !< Particle Jacobian qq
     real(c_real)    :: qp(9)      !< Particle Jacobian qq
     real(c_real)    :: pq(9)      !< Particle Jacobian qq
     real(c_real)    :: pp(9)      !< Particle Jacobian qq
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type fdm_particle_t

end module
