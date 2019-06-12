
# File particle\_mod.f90

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**particle\_mod.f90**](particle__mod_8f90.md)

[Go to the documentation of this file.](particle__mod_8f90.md) 


````cpp
module particle_mod

  use amrex_fort_module, only: c_real => amrex_particle_real
  use iso_c_binding ,    only: c_int

  implicit none
  private

  public  dm_particle_t
  
  type, bind(C)  :: dm_particle_t
     real(c_real)    :: pos(3)
     real(c_real)    :: mass
     real(c_real)    :: vel(3)
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type dm_particle_t

  public  agn_particle_t
  
  type, bind(C)  :: agn_particle_t
     real(c_real)    :: pos(3)
     real(c_real)    :: mass
     real(c_real)    :: vel(3)
     real(c_real)    :: energy
     real(c_real)    :: mdot
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type agn_particle_t

  public  fdm_particle_t
  
  type, bind(C)  :: fdm_particle_t
     real(c_real)    :: pos(3)
     real(c_real)    :: mass
     real(c_real)    :: vel(3)
     real(c_real)    :: phase
     real(c_real)    :: amp(2)
     real(c_real)    :: width
     real(c_real)    :: qq(9)
     real(c_real)    :: pq(9)
     real(c_real)    :: qp(9)
     real(c_real)    :: pp(9)
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type fdm_particle_t

  public  fdm_particle_wkb_t
  
  type, bind(C)  :: fdm_particle_wkb_t
     real(c_real)    :: pos(3)
     real(c_real)    :: mass
     real(c_real)    :: vel(3)
     real(c_real)    :: phase
     real(c_real)    :: amp(2)
     real(c_real)    :: width
     real(c_real)    :: qq(9)
     real(c_real)    :: pq(9)
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type fdm_particle_wkb_t

end module
````

