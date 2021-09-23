#ifndef HYDROHLLE_H
#define HYDROHLLE_H

#include "MFP_hydro_riemann.H"
#include "MFP_hydro_defs.H"

class HydroHLLE : public HydroRiemannSolver
{
public:
    HydroHLLE();
    HydroHLLE(const int i);

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) const override;

    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

};

#endif // HYDROHLLE_H
