#ifndef MFP_HYDRORIEMANN_H
#define MFP_HYDRORIEMANN_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "MFP_factory.H"
#include "MFP_hydro_defs.H"
#include "MFP_utility.H"

using namespace amrex;

//================================================================================

class HydroRiemannSolver
{
public:
    HydroRiemannSolver();
    virtual ~HydroRiemannSolver();

    int idx;

    virtual void solve(Array<Real,+HydroDef::FluxIdx::NUM> &L,
                       Array<Real,+HydroDef::FluxIdx::NUM> &R,
                       Array<Real,+HydroDef::ConsIdx::NUM> &F,
                       Real* shk) {}

    virtual std::string get_tag() const {return "";}
    virtual bool requires_shock_detector() const {return false;}
};

template <typename D>
std::unique_ptr<HydroRiemannSolver> HydroRiemannSolverBuilder(const sol::table& def)
{
    int global_idx = def["global_idx"];

    if ((def["flux"] == D::tag) && (D::valid_state(global_idx))) {
        return std::unique_ptr<D>(new D(global_idx));
    } else {
        return nullptr;
    }
}

ClassFactory<HydroRiemannSolver> &GetHydroRiemannSolverFactory();

#endif // MFP_FLUX_H
