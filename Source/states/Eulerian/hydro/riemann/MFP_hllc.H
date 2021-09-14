#ifndef HYDROHLLC_H
#define HYDROHLLC_H

#include "MFP_hydro_riemann.H"

class HydroHLLC : public HydroRiemannSolver
{
public:
    HydroHLLC();
    HydroHLLC(const int i);

    virtual void solve(Array<Real,+HydroState::FluxIdx::NUM> &L,
                       Array<Real,+HydroState::FluxIdx::NUM> &R,
                       Array<Real,+HydroState::ConsIdx::NUM> &F,
                       Real* shk) const override;

    static bool valid_state(const int idx);
    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

};

#endif // HYDROHLLC_H
