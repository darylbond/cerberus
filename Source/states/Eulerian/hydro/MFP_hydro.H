#ifndef FLUIDSTATE_H
#define FLUIDSTATE_H
#include "forward.hpp"

class HydroState;

#include "MFP_eulerian.H"
#include "MFP_bc.H"
#include "MFP.H"
#include "MFP_hydro_defs.H"
#include "MFP_reconstruction.H"
#include "MFP_hydro_riemann.H"
#include "MFP_hydro_viscous.H"
#include "MFP_hydro_bc.H"
#include "MFP_hydro_shockdetector.H"


class HydroState : public EulerianState
{
public:
    HydroState();
    HydroState(const sol::table& def);
    ~HydroState();

    static bool registered;
    static std::string tag;

    static Vector<std::string> cons_names;
    static Vector<std::string> prim_names;
    static Array<int,1> flux_vector_idx;
    static Array<int,1> cons_vector_idx;
    static Array<int,1> prim_vector_idx;
    static std::map<std::string, int> bc_names;
    static Vector<set_bc> bc_set;

    bool mass_const, charge_const, gamma_const;
    Array<Real,2> mass, charge, gamma;

    std::unique_ptr<HydroRiemannSolver> flux_solver;
    std::unique_ptr<HydroShockDetector> shock_detector;
    std::unique_ptr<HydroViscous> viscous;

    virtual StateType get_type() const override {return StateType::Hydro;}
    virtual bool is_viscous() const {return viscous != nullptr;}

#ifdef AMREX_USE_EB
    Vector<std::unique_ptr<HydroBoundaryEB>> eb_bcs;
    virtual void set_eb_bc(const sol::table &bc_def) override;
#endif

    virtual int get_num_grow() const override {return num_grow;}
    virtual void set_num_grow(int n) override {num_grow = std::max(n, num_grow);}

    virtual int get_num_cons() const override {return +HydroDef::ConsIdx::NUM;}
    virtual const Vector<std::string>& get_cons_names() const override {return cons_names;}

    virtual Vector<std::string> get_plot_output_names() const override;
    virtual void get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                             #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                             #endif
                                 ) const override;

    void set_viscosity();

    Real init_from_number_density(std::map<std::string, Real> data);
    void set_udf();

    void set_flux();
    void set_shock_detector();
    void set_refinement();

    virtual void variable_setup(Vector<int> periodic) override;
    virtual void init_data(MFP* mfp) override;
    virtual void init_from_lua() override;

    Real get_density_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const;
    Real get_density_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q) const;
    Real get_alpha_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const;
    Real get_alpha_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q) const;
    Real get_mass(Real alpha) const;
    Real get_mass(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const;
    Real get_charge(Real alpha) const;
    Real get_charge(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const;
    Real get_gamma(Real alpha) const;
    Real get_gamma(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const;
    Real get_cp(Real alpha) const;
    Real get_cp(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const;
    bool cons2prim(const Array<Real, +HydroDef::ConsIdx::NUM> &U, Array<Real, +HydroDef::PrimIdx::NUM> &Q) const;
    void prim2cons(const Array<Real, +HydroDef::PrimIdx::NUM> &Q, Array<Real, +HydroDef::ConsIdx::NUM> &U) const;
    bool prim_valid(Array<Real,+HydroDef::PrimIdx::NUM>& Q) const;
    bool cons_valid(Array<Real,+HydroDef::ConsIdx::NUM>& U) const;
    Real get_energy_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const;
    Real get_temperature_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const;
    Real get_temperature_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q) const;
    RealArray get_speed_from_cons(const Array<Real, +HydroDef::ConsIdx::NUM> &U) const;
    RealArray get_speed_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q) const;

    virtual Real get_allowed_time_step(MFP* mfp) const override;

    virtual void calc_primitives(const Box& box,
                                 FArrayBox& cons,
                                 FArrayBox &prim,
                                 const Real* dx,
                                 const Real t,
                                 const Real* prob_lo
                             #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                             #endif
                                 ) const override;

    virtual void update_boundary_cells(const Box& box,
                                       const Geometry &geom,
                                       FArrayBox &prim,
                                   #ifdef AMREX_USE_EB
                                       const FArrayBox& vfrac,
                                   #endif
                                       const Real time) const override;

    virtual void calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                 #ifdef AMREX_USE_EB
                                     ,const EBCellFlagFab &flag
                                     ,const FArrayBox &vfrac
                                 #endif
                                     ) const override;

    virtual void calc_time_averaged_faces(const Box& box,
                                          const FArrayBox &prim,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                      #ifdef AMREX_USE_EB
                                          const EBCellFlagFab& flag,
                                      #endif
                                          const Real* dx,
                                          Real dt) const override;

    virtual void face_bc(const int dir,
                         Box const& box,
                         const FArrayBox& src,
                         FArrayBox& dest,
                         const Geometry &geom,
                     #ifdef AMREX_USE_EB
                         const EBCellFlagFab& flag,
                     #endif
                         const Real time,
                         const bool do_all) const override;

    void load_state_for_flux(const Array4<const Real> &face,
                             int i, int j, int k, Array<Real, +HydroDef::FluxIdx::NUM> &S) const;

    virtual void calc_fluxes(const Box& box,
                             FArrayBox &cons,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                             Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                         #ifdef AMREX_USE_EB
                             const EBCellFlagFab& flag,
                         #endif
                             const Real *dx,
                             const Real dt) const override;

    virtual void correct_face_prim(const Box& box,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                                   const Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                               #ifdef AMREX_USE_EB
                                   const EBCellFlagFab &flag,
                               #endif
                                   const Real *dx,
                                   const Real dt) const override;

    void calc_diffusion_terms(const Box& box,
                              const FArrayBox& prim,
                              FArrayBox& diff
                          #ifdef AMREX_USE_EB
                              ,const EBCellFlagFab& flag
                          #endif
                              ) const;

    virtual void calc_viscous_fluxes(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box& pbox,
                                     const FArrayBox &prim,
                                 #ifdef AMREX_USE_EB
                                     const EBCellFlagFab& flag,
                                 #endif
                                     const Real* dx) const override;

#ifdef AMREX_USE_EB

    virtual void calc_viscous_fluxes_eb(const Box& box, Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                        const Box& pbox,
                                        const FArrayBox &prim,
                                    #ifdef AMREX_USE_EB
                                        const EBCellFlagFab& flag,
                                    #endif
                                        const Real* dx) const;

    virtual void calc_wall_fluxes(const Box& box,
                                  const FArrayBox &prim,
                                  Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                  const EBCellFlagFab& flag,
                                  const CutFab &bc_idx,
                                  const FArrayBox& bcent,
                                  const FArrayBox &bnorm,
                                  const Array<const FArrayBox*, AMREX_SPACEDIM> &afrac,
                                  const Real *dx,
                                  const Real dt) const override;
#endif

    virtual void write_info(nlohmann::json& js) const override;

    static HydroState& get_state_global(const int idx) {return static_cast<HydroState&>(MFP::get_state(idx));}

};

#endif // FLUIDSTATE_H
