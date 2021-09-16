#ifndef HYDROHLLE_H
#define HYDROHLLE_H

#include "MFP_hydro_riemann.H"
#include "MFP_hydro_defs.H"

class HydroHLLE : public HydroRiemannSolver
{
public:
    HydroHLLE();
    HydroHLLE(const int i);

    virtual void solve(Array<Real,+HydroDef::FluxIdx::NUM> &L,
                       Array<Real,+HydroDef::FluxIdx::NUM> &R,
                       Array<Real,+HydroDef::ConsIdx::NUM> &F,
                       Real* shk) const override;

    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

};

#endif // HYDROHLLE_H
