#ifndef HYDROHYBRIDHLL_H
#define HYDROHYBRIDHLL_H

#include "MFP_hydro_riemann.H"
#include "MFP_hydro_hlle.H"
#include "MFP_hllc.H"
#include "MFP_hydro_defs.H"


class HydroHybridHLL : public HydroRiemannSolver
{
public:
    HydroHybridHLL();
    HydroHybridHLL(const int i);

    HydroHLLC hllc;
    HydroHLLE hlle;

    virtual void solve(Array<Real,+HydroDef::FluxIdx::NUM> &L,
                       Array<Real,+HydroDef::FluxIdx::NUM> &R,
                       Array<Real,+HydroDef::ConsIdx::NUM> &F,
                       Real* shk) const override;

    virtual bool requires_shock_detector() const {return true;}

    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

private:
    Vector<Real> F_hlle, F_hllc;
};

#endif // HYDROHYBRIDHLL_H
