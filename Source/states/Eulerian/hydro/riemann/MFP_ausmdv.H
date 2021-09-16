#ifndef HYDROAUSMDV_H
#define HYDROAUSMDV_H

#include "MFP_hydro_riemann.H"

class HydroAUSMDV : public HydroRiemannSolver
{
public:
    HydroAUSMDV();
    HydroAUSMDV(const int i);

    virtual void solve(Array<Real,+HydroDef::FluxIdx::NUM> &L,
                       Array<Real,+HydroDef::FluxIdx::NUM> &R,
                       Array<Real,+HydroDef::ConsIdx::NUM> &F,
                       Real* shk) const override;

    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;
};

#endif // HYDROAUSMDV_H
