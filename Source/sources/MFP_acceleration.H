#ifndef ACCELERATION_H
#define ACCELERATION_H

#include "MFP_source.H"
#include "MFP_hydro.H"

class Acceleration : public Source
{
public:
    Acceleration();
    Acceleration(const int idx, const sol::table &def);
    ~Acceleration();

    virtual void solve(MFP* mfp, const Real dt) const override;

    virtual std::string get_tag() const override {return tag;}

    Vector<HydroState*> species;

    Array<Real,3> acc;

    static bool registered;
    static std::string tag;
};

#endif // ACCELERATION_H
