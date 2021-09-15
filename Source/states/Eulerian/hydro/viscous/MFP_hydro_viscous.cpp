#include "MFP_hydro_viscous.H"

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_hydro.H"

// ====================================================================================

HydroViscous::HydroViscous(){}

int HydroViscous::get_type(){return -1;}
int HydroViscous::get_num(){return 0;}


void HydroViscous::update_linked_states()
{
    State& istate = MFP::get_state(idx);
    istate.set_num_grow(2);
}

std::string HydroViscous::str() const
{
    std::stringstream msg;

    const auto coeffs = get_refs();

    msg << get_tag() << "(";

    for (const auto& cf : coeffs) {
        msg << cf.first << "=" << cf.second << ", ";
    }

    msg.seekp(-2,msg.cur);
    msg << ")";

    return msg.str();
}

ClassFactory<HydroViscous>& GetViscousFactory() {
    static ClassFactory<HydroViscous> F;
    return F;
}

// ====================================================================================

Sutherland::Sutherland(){}
Sutherland::~Sutherland(){}

std::string Sutherland::tag = "Sutherland";
bool Sutherland::registered = GetViscousFactory().Register(Sutherland::tag, ViscousBuilder<Sutherland>);

Sutherland::Sutherland(const int global_idx, const sol::table& def)
{

    // check valid state
    if (!valid_state(global_idx))
        Abort("Sutherland viscosity is not valid for "+MFP::get_state(global_idx).name);

    Real mu_ref = def["mu0"];
    Real T_ref = def["T0"];
    Real S_ref = def["S"];
    Prandtl = def["Pr"];
    cfl = def.get_or("cfl",1.0);

    idx = global_idx;

    Real r_r = MFP::rho_ref;
    Real u_r = MFP::u_ref;
    Real x_r = MFP::x_ref;
    Real T_r = MFP::T_ref;

    mu_0 = mu_ref/(r_r*x_r*u_r);
    T0 = T_ref/T_r;
    S = S_ref/T_r;

    if ((mu_0 <= 0) || (T0 <= 0) || (S <= 0)) {
        amrex::Abort("Sutherland input coefficients non-physical");
    }
}

int Sutherland::get_type(){return +HydroViscous::DiffusionType::Neutral;}
int Sutherland::get_num(){return 3;}

void Sutherland::get_coeffs(const Array<Real, +HydroDef::PrimIdx::NUM> &Q, Real &T, Real &mu, Real &kappa)
{
    BL_PROFILE("Sutherland::get_neutral_coeffs");
    const HydroState& istate = HydroState::get_state(idx);

    T = istate.get_temperature_from_prim(Q);
    Real alpha = istate.get_alpha_from_prim(Q);
    Real cp = istate.get_cp(alpha);

    Real T_ = T/T0;
    mu = mu_0*T_*sqrt(T_)*(T0+S)/(T+S);
    kappa = mu*cp/Prandtl;

    return;
}

Real Sutherland::calc_dt(Array<Real,+HydroDef::ConsIdx::NUM>& U, const Real dx2) const
{
    const Real rho = U[+HydroDef::ConsIdx::Density];
    const Real T = HydroState::get_temperature_from_cons(U);
    const Real gamma = HydroState::get_gamma(U);


    const Real T_ = T/T0;
    const Real mu = mu_0*T_*sqrt(T_)*(T0+S)/(T+S);

    return dx2*(4*mu*gamma/(Prandtl*rho));
}

bool Sutherland::valid_state(const int idx)
{

    if (MFP::get_state(idx).get_type() != State::StateType::Hydro) {
        return false;
    }
    return true;
}

// ====================================================================================

PowerLaw::PowerLaw(){}
PowerLaw::~PowerLaw(){}

std::string PowerLaw::tag = "PowerLaw";
bool PowerLaw::registered = GetViscousFactory().Register(PowerLaw::tag, ViscousBuilder<PowerLaw>);

int PowerLaw::get_type(){return +HydroViscous::DiffusionType::Neutral;}
int PowerLaw::get_num(){return 3;}

PowerLaw::PowerLaw(const int global_idx, const sol::table& def)
{

    // check valid state
    if (!valid_state(global_idx))
        Abort("Power Law viscosity is not valid for "+MFP::get_state(global_idx).name);

    Real mu_ref = def["mu0"];
    Real T_ref = def["T0"];
    n = def["n"];
    Prandtl = def["Pr"];
    cfl = def.get_or("cfl",1.0);

    idx = global_idx;

    Real r_r = MFP::rho_ref;
    Real u_r = MFP::u_ref;
    Real x_r = MFP::x_ref;
    Real T_r = MFP::T_ref;

    mu_0 = mu_ref/(r_r*x_r*u_r);
    T0 = T_ref/T_r;


    if ((mu_0 <= 0) || (T0 <= 0) || (n <= 0)) {
        amrex::Abort("Power Law input coefficients non-physical");
    }
}

void PowerLaw::get_coeffs(const Array<Real,+HydroDef::PrimIdx::NUM>& Q, Real &T, Real &mu, Real &kappa)
{
    BL_PROFILE("PowerLaw::get_neutral_coeffs");
    const HydroState& istate = HydroState::get_state(idx);

    T = istate.get_temperature_from_prim(Q);
    Real alpha = istate.get_alpha_from_prim(Q);
    Real cp = istate.get_cp(alpha);


    mu = mu_0*pow(T/T0,n);
    kappa = mu*cp/Prandtl;

    return;
}

Real PowerLaw::calc_dt(Array<Real,+HydroDef::ConsIdx::NUM>& U, const Real dx2) const
{
    const Real rho = U[+HydroDef::ConsIdx::Density];
    const Real T = HydroState::get_temperature_from_cons(U);
    const Real gamma = HydroState::get_gamma(U);

    const Real mu = mu_0*pow(T/T0,n);

    return dx2*(4*mu*gamma/(Prandtl*rho));
}


bool PowerLaw::valid_state(const int idx)
{
    if (MFP::get_state(idx).get_type() != State::StateType::Hydro) {
        return false;
    }
    return true;
}

// ====================================================================================

UserDefinedViscosity::UserDefinedViscosity(){}
UserDefinedViscosity::~UserDefinedViscosity(){}

std::string UserDefinedViscosity::tag = "UserDefined";
bool UserDefinedViscosity::registered = GetViscousFactory().Register(UserDefinedViscosity::tag, ViscousBuilder<UserDefinedViscosity>);

int UserDefinedViscosity::get_type(){return +HydroViscous::DiffusionType::Neutral;}
int UserDefinedViscosity::get_num(){return 3;}

UserDefinedViscosity::UserDefinedViscosity(const int global_idx, const sol::table& def)
{

    // check valid state
    if (!valid_state(global_idx))
        Abort("Constant viscosity is not valid for "+MFP::get_state(global_idx).name);

    mu_0 = def["mu0"];
    Prandtl = def["Pr"];
    cfl = def.get_or("cfl",1.0);

    idx = global_idx;

    if (mu_0 <= 0) {
        amrex::Abort("Constant viscosity input coefficients non-physical");
    }
}

void UserDefinedViscosity::get_coeffs(const Array<Real,+HydroDef::PrimIdx::NUM>& Q, Real &T, Real &mu, Real &kappa)
{
    BL_PROFILE("UserDefinedViscosity::get_neutral_coeffs");
    const HydroState& istate = HydroState::get_state(idx);

    T = istate.get_temperature_from_prim(Q);
    Real alpha = istate.get_alpha_from_prim(Q);
    Real cp = istate.get_cp(alpha);

    mu = mu_0;
    kappa = mu*cp/Prandtl;

    return;
}

Real UserDefinedViscosity::calc_dt(Array<Real,+HydroDef::ConsIdx::NUM>& U, const Real dx2) const
{
    const Real rho = U[+HydroDef::ConsIdx::Density];
    const Real gamma = HydroState::get_gamma(U);
    const Real mu = mu_0;

    return dx2*(4*mu*gamma/(Prandtl*rho));
}

bool UserDefinedViscosity::valid_state(const int idx)
{

    if (MFP::get_state(idx).get_type() != State::StateType::Hydro) {
        return false;
    }
    return true;
}


