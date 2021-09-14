#include "MFP_eulerian.H"
#include "MFP.H"

EulerianState::EulerianState(){}

EulerianState::~EulerianState(){}

void EulerianState::set_reconstruction()
{

    sol::table state_def = MFP::lua["states"][name];
    std::string rec = state_def["reconstruction"].get_or<std::string>("null");

    num_grow = 0;

    return;

}

void EulerianState::set_reflux()
{
    if (AMREX_SPACEDIM < 2) {
        reflux = 0;
    } else {
        reflux = MFP::lua["states"][name]["reflux"].get_or(1);
    }
}

void EulerianState::set_shock_detector()
{


    sol::table sd_def = MFP::lua["states"][name]["shock_detector"].get_or(sol::table());

}


#ifdef AMREX_USE_EB
void EulerianState::set_eb_divergence()
{

    sol::table state_def = MFP::lua["states"][name];


    sol::optional<sol::table> div_def_exists = state_def["eb_divergence"];

    // set a default
    if (!div_def_exists) {
        state_def["eb_divergence"] = MFP::lua["eb_divergence_default"];
    }

    sol::table div_def = state_def["eb_divergence"];

    return;
}

#endif

void EulerianState::init_from_lua()
{
    BL_PROFILE("EulerianState::init_from_lua");

    //
    // reflux
    //
    set_reflux();

    //
    // get reconstruction
    //

    set_reconstruction();

#ifdef AMREX_USE_EB
    //
    // how to handle boundary divergence
    //
    set_eb_divergence();
#endif

}

