#include "MFP_eulerian.H"
#include "MFP.H"
#include "MFP_shockdetector.H"

EulerianState::EulerianState(){}

EulerianState::~EulerianState(){}

void EulerianState::set_reconstruction()
{
    ClassFactory<Reconstruction> rfact = GetReconstructionFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string rec = state_def["reconstruction"].get_or<std::string>("null");

    state_def["reconstruction"] = rec; // consistency when using default option "null"

    if (is_transported() && rec == "null")
        Abort("Reconstruction option required for state '"+name+"'. Options are "+vec2str(rfact.getKeys()));

    reconstruction = rfact.Build(rec, state_def);

    if (!reconstruction)
        Abort("Invalid reconstruction option '"+rec+"'. Options are "+vec2str(rfact.getKeys()));

    set_num_grow(reconstruction->get_num_grow());

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

    ClassFactory<ShockDetector> sdfact = GetShockDetectorFactory();

    sol::table sd_def = MFP::lua["states"][name]["shock_detector"].get_or(sol::table());

    if (!sd_def.valid()) return;

    sd_def["global_idx"] = global_idx;

    std::string sd_name = sd_def["name"].get_or<std::string>("");

    shock_detector = sdfact.Build(sd_name, sd_def);

    if (!sd_name.empty() && !shock_detector)
        Abort("Invalid shock_detector option '"+sd_name+"'. Options are "+vec2str(sdfact.getKeys()));
}


#ifdef AMREX_USE_EB
void EulerianState::set_eb_divergence()
{

    ClassFactory<DivergenceEB> div_fact = GetDivergenceEBBuilder();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;

    sol::optional<sol::table> div_def_exists = state_def["eb_divergence"];

    // set a default
    if (!div_def_exists) {
        state_def["eb_divergence"] = MFP::lua["eb_divergence_default"];
    }

    sol::table div_def = state_def["eb_divergence"];

    std::string tag = div_def.get<std::string>("type");
    eb_div = div_fact.Build(tag, state_def);

    if (!eb_div)
        Abort("Invalid 'eb_divergence' type '"+tag+"'. Options are "+vec2str(div_fact.getKeys()));

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

