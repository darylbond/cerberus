#include "MFP_shockdetector.H"

ShockDetector::ShockDetector(){}

ShockDetector::~ShockDetector(){}

Real ShockDetector::solve(Vector<Real> &L,Vector<Real> &R) const {return 0.0;}

ClassFactory<ShockDetector>& GetShockDetectorFactory()
{
    static ClassFactory<ShockDetector> F;
    return F;
}

//================================================================================

#include "MFP_hydro.H"

std::string HydroShockDetector::tag = "pressure_jump_detector";
bool HydroShockDetector::registered = GetShockDetectorFactory().Register(HydroShockDetector::tag, ShockDetectorBuilder<HydroShockDetector>);

HydroShockDetector::HydroShockDetector(){}
HydroShockDetector::HydroShockDetector(const sol::table& def)
{
    idx = def["global_idx"];

    shock_threshold = def["threshold"].get_or(-1.0);

    if (shock_threshold < 0.0) {
        Abort("Error: 'pressure_jump_detector' requires 'threshold' to be set (>0.0)");
    }

}

Real HydroShockDetector::solve(Vector<Real> &L, Vector<Real> &R) const
{
    BL_PROFILE("HydroShockDetector::solve");

    Real pL = L[+HydroDef::PrimIdx::Prs];
    Real pR = R[+HydroDef::PrimIdx::Prs];
    Real varphi = std::abs(pR - pL)/(pR + pL);

    return 0.5 + 0.5*tanh_approx(5*(varphi - 0.75*shock_threshold)/shock_threshold);
}

bool HydroShockDetector::valid_state(const int idx)
{

    if (MFP::get_state(idx).get_type() == State::StateType::Hydro ) return true;

    return false;
}

void HydroShockDetector::write_info(nlohmann::json& js) const
{
    nlohmann::json& sd = js["shock_detector"];

    sd["name"] = tag;
    sd["threshold"] = shock_threshold;
}
