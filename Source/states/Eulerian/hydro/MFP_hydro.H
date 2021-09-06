#ifndef FLUIDSTATE_H
#define FLUIDSTATE_H
#include "forward.hpp"

#include "MFP_eulerian.H"
#include "MFP_bc.H"
#include "MFP.H"

class HydroRiemannSolver;
#include "MFP_hydro_riemann.H"

class HydroState;
class HydroViscous;

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
    int enforce_positivity;
    int extra_slope_limits;

    std::unique_ptr<HydroViscous> viscous;
    std::unique_ptr<HydroRiemannSolver> flux_solver;
    std::unique_ptr<ShockDetector> shock_detector;

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
#ifdef AMREX_USE_EB
    virtual void set_eb_bc(const sol::table &bc_def) override;
#endif
    virtual void init_from_lua() override;

    static Real get_density_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U);
    static Real get_density_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q);
    static Real get_alpha_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U);
    static Real get_alpha_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q);
    static Real get_mass(Real alpha);
    static Real get_mass(const Array<Real,+HydroDef::ConsIdx::NUM>& U);
    static Real get_charge(Real alpha);
    static Real get_charge(const Array<Real,+HydroDef::ConsIdx::NUM>& U);
    static Real get_gamma(Real alpha);
    static Real get_gamma(const Array<Real,+HydroDef::ConsIdx::NUM>& U);
    static Real get_cp(Real alpha);
    static Real get_cp(const Array<Real,+HydroDef::ConsIdx::NUM>& U);
    static bool cons2prim(Vector<Real>& U, Vector<Real>& Q);
    static void prim2cons(Vector<Real>& Q, Vector<Real>& U);
    static bool prim_valid(Array<Real,+HydroDef::PrimIdx::NUM>& Q);
    static bool cons_valid(Array<Real,+HydroDef::ConsIdx::NUM>& U);
    static Real get_energy_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U);
    static Real get_temperature_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U);
    static Real get_temperature_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q);
    static RealArray get_speed_from_cons(const Array<Real, +HydroDef::ConsIdx::NUM> &U);
    static RealArray get_speed_from_prim(const Vector<Real>& Q);

    virtual Real get_allowed_time_step(MFP* mfp) const override;


    static HydroState& get_state(const int idx) {return static_cast<HydroState&>(MFP::get_state(idx));}

};

#endif // FLUIDSTATE_H
