
# Class Gravity


[**Class List**](annotated.md) **>** [**Gravity**](classGravity.md)


















## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::array&lt; amrex::MLLinOp::BCType, AMREX\_SPACEDIM &gt; | [**mlmg\_hibc**](classGravity.md#variable-mlmg-hibc)  <br> |
|  std::array&lt; amrex::MLLinOp::BCType, AMREX\_SPACEDIM &gt; | [**mlmg\_lobc**](classGravity.md#variable-mlmg-lobc)  <br> |


## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Gravity**](classGravity.md#function-gravity) (amrex::Amr \* Parent, int \_finest\_level, amrex::BCRec \* \_phys\_bc, int \_Density) <br> |
|  void | [**actual\_multilevel\_solve**](classGravity.md#function-actual-multilevel-solve) (int level, int finest\_level, const amrex::Vector&lt; amrex::Vector&lt; amrex::MultiFab \* &gt; &gt; & grad\_phi, int is\_new, int ngrow\_for\_solve, int use\_previous\_phi\_as\_guess=0) <br> |
|  void | [**add\_to\_fluxes**](classGravity.md#function-add-to-fluxes) (int level, int iteration, int ncycle) <br> |
|  void | [**average\_fine\_ec\_onto\_crse\_ec**](classGravity.md#function-average-fine-ec-onto-crse-ec) (int level, int is\_new) <br> |
|  amrex::Real | [**compute\_level\_average**](classGravity.md#function-compute-level-average) (int level, amrex::MultiFab \* mf) <br> |
|  amrex::Real | [**compute\_multilevel\_average**](classGravity.md#function-compute-multilevel-average) (int level, amrex::MultiFab \* mf, int flev=-1) <br> |
|  amrex::Real | [**get\_const\_grav**](classGravity.md#function-get-const-grav) () <br> |
|  void | [**get\_crse\_grad\_phi**](classGravity.md#function-get-crse-grad-phi) (int level, amrex::Vector&lt; std::unique\_ptr&lt; amrex::MultiFab &gt; &gt; & grad\_phi\_crse, amrex::Real time) <br> |
|  void | [**get\_crse\_phi**](classGravity.md#function-get-crse-phi) (int level, amrex::MultiFab & phi\_crse, amrex::Real time) <br> |
|  amrex::Vector&lt; amrex::MultiFab \* &gt; | [**get\_grad\_phi\_curr**](classGravity.md#function-get-grad-phi-curr) (int level) <br> |
|  amrex::Vector&lt; amrex::MultiFab \* &gt; | [**get\_grad\_phi\_prev**](classGravity.md#function-get-grad-phi-prev) (int level) <br> |
|  std::string | [**get\_gravity\_type**](classGravity.md#function-get-gravity-type) () <br> |
|  void | [**get\_new\_grav\_vector**](classGravity.md#function-get-new-grav-vector) (int level, amrex::MultiFab & grav\_vector, amrex::Real time) <br> |
|  int | [**get\_no\_composite**](classGravity.md#function-get-no-composite) () <br> |
|  int | [**get\_no\_sync**](classGravity.md#function-get-no-sync) () <br> |
|  void | [**get\_old\_grav\_vector**](classGravity.md#function-get-old-grav-vector) (int level, amrex::MultiFab & grav\_vector, amrex::Real time) <br> |
|  void | [**gravity\_sync**](classGravity.md#function-gravity-sync) (int crse\_level, int fine\_level, int iteration, int ncycle, const amrex::MultiFab & drho\_and\_drhoU, const amrex::MultiFab & dphi, const amrex::Vector&lt; amrex::MultiFab \* &gt; & grad\_delta\_phi\_cc) <br> |
|  void | [**install\_level**](classGravity.md#function-install-level) (int level, amrex::AmrLevel \* level\_data\_to\_install) <br> |
|  void | [**make\_mg\_bc**](classGravity.md#function-make-mg-bc) () <br> |
|  void | [**multilevel\_solve\_for\_new\_phi**](classGravity.md#function-multilevel-solve-for-new-phi) (int level, int finest\_level, int ngrow\_for\_solve, int use\_previous\_phi\_as\_guess=0) <br> |
|  void | [**multilevel\_solve\_for\_old\_phi**](classGravity.md#function-multilevel-solve-for-old-phi) (int level, int finest\_level, int ngrow\_for\_solve, int use\_previous\_phi\_as\_guess=0) <br> |
|  void | [**plus\_grad\_phi\_curr**](classGravity.md#function-plus-grad-phi-curr) (int level, const amrex::Vector&lt; amrex::MultiFab \* &gt; & addend) <br> |
|  void | [**plus\_phi\_curr**](classGravity.md#function-plus-phi-curr) (int level, amrex::MultiFab & addend) <br> |
|  void | [**read\_params**](classGravity.md#function-read-params) () <br> |
|  void | [**reflux\_phi**](classGravity.md#function-reflux-phi) (int level, amrex::MultiFab & dphi) <br> |
|  void | [**set\_boundary**](classGravity.md#function-set-boundary) (amrex::BndryData & bd, amrex::MultiFab & rhs, const amrex::Real \* dx) <br> |
|  void | [**set\_dirichlet\_bcs**](classGravity.md#function-set-dirichlet-bcs) (int level, amrex::MultiFab \* phi) <br> |
|  void | [**set\_mass\_offset**](classGravity.md#function-set-mass-offset) (amrex::Real time) <br> |
|  void | [**solve\_for\_delta\_phi**](classGravity.md#function-solve-for-delta-phi) (int crse\_level, int fine\_level, amrex::MultiFab & CrseRhs, const amrex::Vector&lt; amrex::MultiFab \* &gt; & delta\_phi, const amrex::Vector&lt; amrex::Vector&lt; amrex::MultiFab \* &gt; &gt; & grad\_delta\_phi) <br> |
|  void | [**solve\_for\_new\_phi**](classGravity.md#function-solve-for-new-phi) (int level, amrex::MultiFab & phi, const amrex::Vector&lt; amrex::MultiFab \* &gt; & grad\_phi, int fill\_interior, int grav\_n\_grow=1) <br> |
|  void | [**solve\_for\_old\_phi**](classGravity.md#function-solve-for-old-phi) (int level, amrex::MultiFab & phi, const amrex::Vector&lt; amrex::MultiFab \* &gt; & grad\_phi, int fill\_interior, int grav\_n\_grow=1) <br> |
|  void | [**solve\_for\_phi**](classGravity.md#function-solve-for-phi) (int level, amrex::MultiFab & Rhs, amrex::MultiFab & phi, const amrex::Vector&lt; amrex::MultiFab \* &gt; & grad\_phi, amrex::Real time, int fill\_interior) <br> |
|  amrex::Real | [**solve\_with\_MLMG**](classGravity.md#function-solve-with-mlmg) (int crse\_level, int fine\_level, const amrex::Vector&lt; amrex::MultiFab \* &gt; & phi, const amrex::Vector&lt; const amrex::MultiFab \* &gt; & rhs, const amrex::Vector&lt; std::array&lt; amrex::MultiFab \*, AMREX\_SPACEDIM &gt; &gt; & grad\_phi, const amrex::MultiFab \*const crse\_bcdata, amrex::Real rel\_eps, amrex::Real abs\_eps) <br> |
|  void | [**swap\_time\_levels**](classGravity.md#function-swap-time-levels) (int level) <br> |
|  void | [**zero\_phi\_flux\_reg**](classGravity.md#function-zero-phi-flux-reg) (int level) <br> |
| virtual  | [**~Gravity**](classGravity.md#function-gravity) () <br> |




## Protected Attributes

| Type | Name |
| ---: | :--- |
|  amrex::Vector&lt; amrex::AmrLevel \* &gt; | [**LevelData**](classGravity.md#variable-leveldata)  <br> |
|  int | [**density**](classGravity.md#variable-density)  <br> |
|  const amrex::Vector&lt; amrex::DistributionMapping &gt; & | [**dmap**](classGravity.md#variable-dmap)  <br> |
|  int | [**finest\_level**](classGravity.md#variable-finest-level)  <br> |
|  int | [**finest\_level\_allocated**](classGravity.md#variable-finest-level-allocated)  <br> |
|  amrex::Vector&lt; amrex::Vector&lt; std::unique\_ptr&lt; amrex::MultiFab &gt; &gt; &gt; | [**grad\_phi\_curr**](classGravity.md#variable-grad-phi-curr)  <br> |
|  amrex::Vector&lt; amrex::Vector&lt; std::unique\_ptr&lt; amrex::MultiFab &gt; &gt; &gt; | [**grad\_phi\_prev**](classGravity.md#variable-grad-phi-prev)  <br> |
|  const amrex::Vector&lt; amrex::BoxArray &gt; & | [**grids**](classGravity.md#variable-grids)  <br> |
|  amrex::Vector&lt; amrex::Real &gt; | [**level\_solver\_resnorm**](classGravity.md#variable-level-solver-resnorm)  <br> |
|  amrex::Amr \* | [**parent**](classGravity.md#variable-parent)  <br> |
|  amrex::Vector&lt; std::unique\_ptr&lt; amrex::FluxRegister &gt; &gt; | [**phi\_flux\_reg**](classGravity.md#variable-phi-flux-reg)  <br> |
|  amrex::BCRec \* | [**phys\_bc**](classGravity.md#variable-phys-bc)  <br> |

## Protected Static Attributes

| Type | Name |
| ---: | :--- |
|  amrex::Real | [**delta\_tol**](classGravity.md#variable-delta-tol)   = = 1.e-12<br> |
|  int | [**dirichlet\_bcs**](classGravity.md#variable-dirichlet-bcs)   = = 0<br> |
|  std::string | [**gravity\_type**](classGravity.md#variable-gravity-type)   = = "fill\_me\_please"<br> |
|  amrex::Real | [**mass\_offset**](classGravity.md#variable-mass-offset)   = = 0<br> |
|  amrex::Real | [**ml\_tol**](classGravity.md#variable-ml-tol)   = = 1.e-12<br> |
|  int | [**mlmg\_agglomeration**](classGravity.md#variable-mlmg-agglomeration)   = = 0<br> |
|  int | [**mlmg\_consolidation**](classGravity.md#variable-mlmg-consolidation)   = = 0<br> |
|  int | [**mlmg\_max\_fmg\_iter**](classGravity.md#variable-mlmg-max-fmg-iter)   = = 0<br> |
|  int | [**no\_composite**](classGravity.md#variable-no-composite)   = = 0<br> |
|  int | [**no\_sync**](classGravity.md#variable-no-sync)   = = 0<br> |
|  amrex::Real | [**sl\_tol**](classGravity.md#variable-sl-tol)   = = 1.e-12<br> |
|  int | [**verbose**](classGravity.md#variable-verbose)   = = 0<br> |

## Protected Functions

| Type | Name |
| ---: | :--- |
|  void | [**AddGhostParticlesToRhs**](classGravity.md#function-addghostparticlestorhs-1-2) (int level, amrex::MultiFab & Rhs, int ngrow) <br> |
|  void | [**AddGhostParticlesToRhs**](classGravity.md#function-addghostparticlestorhs-2-2) (int level, const amrex::Vector&lt; amrex::MultiFab \* &gt; & Rhs\_particles, int ngrow) <br> |
|  void | [**AddParticlesToRhs**](classGravity.md#function-addparticlestorhs-1-2) (int level, amrex::MultiFab & Rhs, int ngrow) <br> |
|  void | [**AddParticlesToRhs**](classGravity.md#function-addparticlestorhs-2-2) (int base\_level, int finest\_level, int ngrow, const amrex::Vector&lt; amrex::MultiFab \* &gt; & Rhs\_particles) <br> |
|  void | [**AddVirtualParticlesToRhs**](classGravity.md#function-addvirtualparticlestorhs-1-2) (int level, amrex::MultiFab & Rhs, int ngrow) <br> |
|  void | [**AddVirtualParticlesToRhs**](classGravity.md#function-addvirtualparticlestorhs-2-2) (int finest\_level, const amrex::Vector&lt; amrex::MultiFab \* &gt; & Rhs\_particles, int ngrow) <br> |
|  void | [**CorrectRhsUsingOffset**](classGravity.md#function-correctrhsusingoffset) (int level, amrex::MultiFab & Rhs) <br> |
|  void | [**fill\_ec\_grow**](classGravity.md#function-fill-ec-grow) (int level, const amrex::Vector&lt; amrex::MultiFab \* &gt; & ecF, const amrex::Vector&lt; amrex::MultiFab \* &gt; & ecC) const<br> |


## Public Attributes Documentation


### <a href="#variable-mlmg-hibc" id="variable-mlmg-hibc">variable mlmg\_hibc </a>


```cpp
std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> Gravity::mlmg_hibc;
```



### <a href="#variable-mlmg-lobc" id="variable-mlmg-lobc">variable mlmg\_lobc </a>


```cpp
std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> Gravity::mlmg_lobc;
```


## Public Functions Documentation


### <a href="#function-gravity" id="function-gravity">function Gravity </a>


```cpp
Gravity::Gravity (
    amrex::Amr * Parent,
    int _finest_level,
    amrex::BCRec * _phys_bc,
    int _Density
) 
```



### <a href="#function-actual-multilevel-solve" id="function-actual-multilevel-solve">function actual\_multilevel\_solve </a>


```cpp
void Gravity::actual_multilevel_solve (
    int level,
    int finest_level,
    const amrex::Vector< amrex::Vector< amrex::MultiFab * > > & grad_phi,
    int is_new,
    int ngrow_for_solve,
    int use_previous_phi_as_guess=0
) 
```



### <a href="#function-add-to-fluxes" id="function-add-to-fluxes">function add\_to\_fluxes </a>


```cpp
void Gravity::add_to_fluxes (
    int level,
    int iteration,
    int ncycle
) 
```



### <a href="#function-average-fine-ec-onto-crse-ec" id="function-average-fine-ec-onto-crse-ec">function average\_fine\_ec\_onto\_crse\_ec </a>


```cpp
void Gravity::average_fine_ec_onto_crse_ec (
    int level,
    int is_new
) 
```



### <a href="#function-compute-level-average" id="function-compute-level-average">function compute\_level\_average </a>


```cpp
amrex::Real Gravity::compute_level_average (
    int level,
    amrex::MultiFab * mf
) 
```



### <a href="#function-compute-multilevel-average" id="function-compute-multilevel-average">function compute\_multilevel\_average </a>


```cpp
amrex::Real Gravity::compute_multilevel_average (
    int level,
    amrex::MultiFab * mf,
    int flev=-1
) 
```



### <a href="#function-get-const-grav" id="function-get-const-grav">function get\_const\_grav </a>


```cpp
amrex::Real Gravity::get_const_grav () 
```



### <a href="#function-get-crse-grad-phi" id="function-get-crse-grad-phi">function get\_crse\_grad\_phi </a>


```cpp
void Gravity::get_crse_grad_phi (
    int level,
    amrex::Vector< std::unique_ptr< amrex::MultiFab > > & grad_phi_crse,
    amrex::Real time
) 
```



### <a href="#function-get-crse-phi" id="function-get-crse-phi">function get\_crse\_phi </a>


```cpp
void Gravity::get_crse_phi (
    int level,
    amrex::MultiFab & phi_crse,
    amrex::Real time
) 
```



### <a href="#function-get-grad-phi-curr" id="function-get-grad-phi-curr">function get\_grad\_phi\_curr </a>


```cpp
amrex::Vector< amrex::MultiFab * > Gravity::get_grad_phi_curr (
    int level
) 
```



### <a href="#function-get-grad-phi-prev" id="function-get-grad-phi-prev">function get\_grad\_phi\_prev </a>


```cpp
amrex::Vector< amrex::MultiFab * > Gravity::get_grad_phi_prev (
    int level
) 
```



### <a href="#function-get-gravity-type" id="function-get-gravity-type">function get\_gravity\_type </a>


```cpp
std::string Gravity::get_gravity_type () 
```



### <a href="#function-get-new-grav-vector" id="function-get-new-grav-vector">function get\_new\_grav\_vector </a>


```cpp
void Gravity::get_new_grav_vector (
    int level,
    amrex::MultiFab & grav_vector,
    amrex::Real time
) 
```



### <a href="#function-get-no-composite" id="function-get-no-composite">function get\_no\_composite </a>


```cpp
int Gravity::get_no_composite () 
```



### <a href="#function-get-no-sync" id="function-get-no-sync">function get\_no\_sync </a>


```cpp
int Gravity::get_no_sync () 
```



### <a href="#function-get-old-grav-vector" id="function-get-old-grav-vector">function get\_old\_grav\_vector </a>


```cpp
void Gravity::get_old_grav_vector (
    int level,
    amrex::MultiFab & grav_vector,
    amrex::Real time
) 
```



### <a href="#function-gravity-sync" id="function-gravity-sync">function gravity\_sync </a>


```cpp
void Gravity::gravity_sync (
    int crse_level,
    int fine_level,
    int iteration,
    int ncycle,
    const amrex::MultiFab & drho_and_drhoU,
    const amrex::MultiFab & dphi,
    const amrex::Vector< amrex::MultiFab * > & grad_delta_phi_cc
) 
```



### <a href="#function-install-level" id="function-install-level">function install\_level </a>


```cpp
void Gravity::install_level (
    int level,
    amrex::AmrLevel * level_data_to_install
) 
```



### <a href="#function-make-mg-bc" id="function-make-mg-bc">function make\_mg\_bc </a>


```cpp
void Gravity::make_mg_bc () 
```



### <a href="#function-multilevel-solve-for-new-phi" id="function-multilevel-solve-for-new-phi">function multilevel\_solve\_for\_new\_phi </a>


```cpp
void Gravity::multilevel_solve_for_new_phi (
    int level,
    int finest_level,
    int ngrow_for_solve,
    int use_previous_phi_as_guess=0
) 
```



### <a href="#function-multilevel-solve-for-old-phi" id="function-multilevel-solve-for-old-phi">function multilevel\_solve\_for\_old\_phi </a>


```cpp
void Gravity::multilevel_solve_for_old_phi (
    int level,
    int finest_level,
    int ngrow_for_solve,
    int use_previous_phi_as_guess=0
) 
```



### <a href="#function-plus-grad-phi-curr" id="function-plus-grad-phi-curr">function plus\_grad\_phi\_curr </a>


```cpp
void Gravity::plus_grad_phi_curr (
    int level,
    const amrex::Vector< amrex::MultiFab * > & addend
) 
```



### <a href="#function-plus-phi-curr" id="function-plus-phi-curr">function plus\_phi\_curr </a>


```cpp
void Gravity::plus_phi_curr (
    int level,
    amrex::MultiFab & addend
) 
```



### <a href="#function-read-params" id="function-read-params">function read\_params </a>


```cpp
void Gravity::read_params () 
```



### <a href="#function-reflux-phi" id="function-reflux-phi">function reflux\_phi </a>


```cpp
void Gravity::reflux_phi (
    int level,
    amrex::MultiFab & dphi
) 
```



### <a href="#function-set-boundary" id="function-set-boundary">function set\_boundary </a>


```cpp
void Gravity::set_boundary (
    amrex::BndryData & bd,
    amrex::MultiFab & rhs,
    const amrex::Real * dx
) 
```



### <a href="#function-set-dirichlet-bcs" id="function-set-dirichlet-bcs">function set\_dirichlet\_bcs </a>


```cpp
void Gravity::set_dirichlet_bcs (
    int level,
    amrex::MultiFab * phi
) 
```



### <a href="#function-set-mass-offset" id="function-set-mass-offset">function set\_mass\_offset </a>


```cpp
void Gravity::set_mass_offset (
    amrex::Real time
) 
```



### <a href="#function-solve-for-delta-phi" id="function-solve-for-delta-phi">function solve\_for\_delta\_phi </a>


```cpp
void Gravity::solve_for_delta_phi (
    int crse_level,
    int fine_level,
    amrex::MultiFab & CrseRhs,
    const amrex::Vector< amrex::MultiFab * > & delta_phi,
    const amrex::Vector< amrex::Vector< amrex::MultiFab * > > & grad_delta_phi
) 
```



### <a href="#function-solve-for-new-phi" id="function-solve-for-new-phi">function solve\_for\_new\_phi </a>


```cpp
void Gravity::solve_for_new_phi (
    int level,
    amrex::MultiFab & phi,
    const amrex::Vector< amrex::MultiFab * > & grad_phi,
    int fill_interior,
    int grav_n_grow=1
) 
```



### <a href="#function-solve-for-old-phi" id="function-solve-for-old-phi">function solve\_for\_old\_phi </a>


```cpp
void Gravity::solve_for_old_phi (
    int level,
    amrex::MultiFab & phi,
    const amrex::Vector< amrex::MultiFab * > & grad_phi,
    int fill_interior,
    int grav_n_grow=1
) 
```



### <a href="#function-solve-for-phi" id="function-solve-for-phi">function solve\_for\_phi </a>


```cpp
void Gravity::solve_for_phi (
    int level,
    amrex::MultiFab & Rhs,
    amrex::MultiFab & phi,
    const amrex::Vector< amrex::MultiFab * > & grad_phi,
    amrex::Real time,
    int fill_interior
) 
```



### <a href="#function-solve-with-mlmg" id="function-solve-with-mlmg">function solve\_with\_MLMG </a>


```cpp
amrex::Real Gravity::solve_with_MLMG (
    int crse_level,
    int fine_level,
    const amrex::Vector< amrex::MultiFab * > & phi,
    const amrex::Vector< const amrex::MultiFab * > & rhs,
    const amrex::Vector< std::array< amrex::MultiFab *, AMREX_SPACEDIM > > & grad_phi,
    const amrex::MultiFab *const crse_bcdata,
    amrex::Real rel_eps,
    amrex::Real abs_eps
) 
```



### <a href="#function-swap-time-levels" id="function-swap-time-levels">function swap\_time\_levels </a>


```cpp
void Gravity::swap_time_levels (
    int level
) 
```



### <a href="#function-zero-phi-flux-reg" id="function-zero-phi-flux-reg">function zero\_phi\_flux\_reg </a>


```cpp
void Gravity::zero_phi_flux_reg (
    int level
) 
```



### <a href="#function-gravity" id="function-gravity">function ~Gravity </a>


```cpp
virtual Gravity::~Gravity () 
```


## Protected Attributes Documentation


### <a href="#variable-leveldata" id="variable-leveldata">variable LevelData </a>


```cpp
amrex::Vector<amrex::AmrLevel*> Gravity::LevelData;
```



### <a href="#variable-density" id="variable-density">variable density </a>


```cpp
int Gravity::density;
```



### <a href="#variable-dmap" id="variable-dmap">variable dmap </a>


```cpp
const amrex::Vector<amrex::DistributionMapping>& Gravity::dmap;
```



### <a href="#variable-finest-level" id="variable-finest-level">variable finest\_level </a>


```cpp
int Gravity::finest_level;
```



### <a href="#variable-finest-level-allocated" id="variable-finest-level-allocated">variable finest\_level\_allocated </a>


```cpp
int Gravity::finest_level_allocated;
```



### <a href="#variable-grad-phi-curr" id="variable-grad-phi-curr">variable grad\_phi\_curr </a>


```cpp
amrex::Vector< amrex::Vector<std::unique_ptr<amrex::MultiFab> > > Gravity::grad_phi_curr;
```



### <a href="#variable-grad-phi-prev" id="variable-grad-phi-prev">variable grad\_phi\_prev </a>


```cpp
amrex::Vector< amrex::Vector<std::unique_ptr<amrex::MultiFab> > > Gravity::grad_phi_prev;
```



### <a href="#variable-grids" id="variable-grids">variable grids </a>


```cpp
const amrex::Vector<amrex::BoxArray>& Gravity::grids;
```



### <a href="#variable-level-solver-resnorm" id="variable-level-solver-resnorm">variable level\_solver\_resnorm </a>


```cpp
amrex::Vector<amrex::Real> Gravity::level_solver_resnorm;
```



### <a href="#variable-parent" id="variable-parent">variable parent </a>


```cpp
amrex::Amr* Gravity::parent;
```



### <a href="#variable-phi-flux-reg" id="variable-phi-flux-reg">variable phi\_flux\_reg </a>


```cpp
amrex::Vector<std::unique_ptr<amrex::FluxRegister> > Gravity::phi_flux_reg;
```



### <a href="#variable-phys-bc" id="variable-phys-bc">variable phys\_bc </a>


```cpp
amrex::BCRec* Gravity::phys_bc;
```


## Protected Static Attributes Documentation


### <a href="#variable-delta-tol" id="variable-delta-tol">variable delta\_tol </a>


```cpp
Real Gravity::delta_tol;
```



### <a href="#variable-dirichlet-bcs" id="variable-dirichlet-bcs">variable dirichlet\_bcs </a>


```cpp
int Gravity::dirichlet_bcs;
```



### <a href="#variable-gravity-type" id="variable-gravity-type">variable gravity\_type </a>


```cpp
std::string Gravity::gravity_type;
```



### <a href="#variable-mass-offset" id="variable-mass-offset">variable mass\_offset </a>


```cpp
Real Gravity::mass_offset;
```



### <a href="#variable-ml-tol" id="variable-ml-tol">variable ml\_tol </a>


```cpp
Real Gravity::ml_tol;
```



### <a href="#variable-mlmg-agglomeration" id="variable-mlmg-agglomeration">variable mlmg\_agglomeration </a>


```cpp
int Gravity::mlmg_agglomeration;
```



### <a href="#variable-mlmg-consolidation" id="variable-mlmg-consolidation">variable mlmg\_consolidation </a>


```cpp
int Gravity::mlmg_consolidation;
```



### <a href="#variable-mlmg-max-fmg-iter" id="variable-mlmg-max-fmg-iter">variable mlmg\_max\_fmg\_iter </a>


```cpp
int Gravity::mlmg_max_fmg_iter;
```



### <a href="#variable-no-composite" id="variable-no-composite">variable no\_composite </a>


```cpp
int Gravity::no_composite;
```



### <a href="#variable-no-sync" id="variable-no-sync">variable no\_sync </a>


```cpp
int Gravity::no_sync;
```



### <a href="#variable-sl-tol" id="variable-sl-tol">variable sl\_tol </a>


```cpp
Real Gravity::sl_tol;
```



### <a href="#variable-verbose" id="variable-verbose">variable verbose </a>


```cpp
int Gravity::verbose;
```


## Protected Functions Documentation


### <a href="#function-addghostparticlestorhs-1-2" id="function-addghostparticlestorhs-1-2">function AddGhostParticlesToRhs [1/2]</a>


```cpp
void Gravity::AddGhostParticlesToRhs (
    int level,
    amrex::MultiFab & Rhs,
    int ngrow
) 
```



### <a href="#function-addghostparticlestorhs-2-2" id="function-addghostparticlestorhs-2-2">function AddGhostParticlesToRhs [2/2]</a>


```cpp
void Gravity::AddGhostParticlesToRhs (
    int level,
    const amrex::Vector< amrex::MultiFab * > & Rhs_particles,
    int ngrow
) 
```



### <a href="#function-addparticlestorhs-1-2" id="function-addparticlestorhs-1-2">function AddParticlesToRhs [1/2]</a>


```cpp
void Gravity::AddParticlesToRhs (
    int level,
    amrex::MultiFab & Rhs,
    int ngrow
) 
```



### <a href="#function-addparticlestorhs-2-2" id="function-addparticlestorhs-2-2">function AddParticlesToRhs [2/2]</a>


```cpp
void Gravity::AddParticlesToRhs (
    int base_level,
    int finest_level,
    int ngrow,
    const amrex::Vector< amrex::MultiFab * > & Rhs_particles
) 
```



### <a href="#function-addvirtualparticlestorhs-1-2" id="function-addvirtualparticlestorhs-1-2">function AddVirtualParticlesToRhs [1/2]</a>


```cpp
void Gravity::AddVirtualParticlesToRhs (
    int level,
    amrex::MultiFab & Rhs,
    int ngrow
) 
```



### <a href="#function-addvirtualparticlestorhs-2-2" id="function-addvirtualparticlestorhs-2-2">function AddVirtualParticlesToRhs [2/2]</a>


```cpp
void Gravity::AddVirtualParticlesToRhs (
    int finest_level,
    const amrex::Vector< amrex::MultiFab * > & Rhs_particles,
    int ngrow
) 
```



### <a href="#function-correctrhsusingoffset" id="function-correctrhsusingoffset">function CorrectRhsUsingOffset </a>


```cpp
void Gravity::CorrectRhsUsingOffset (
    int level,
    amrex::MultiFab & Rhs
) 
```



### <a href="#function-fill-ec-grow" id="function-fill-ec-grow">function fill\_ec\_grow </a>


```cpp
void Gravity::fill_ec_grow (
    int level,
    const amrex::Vector< amrex::MultiFab * > & ecF,
    const amrex::Vector< amrex::MultiFab * > & ecC
) const
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Gravity/Gravity.H`