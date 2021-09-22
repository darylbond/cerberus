#ifndef PLASMA5_H
#define PLASMA5_H

#include "MFP_source.H"
#include "MFP_field.H"
#include "MFP_hydro.H"

class Plasma5 : public Source
{
public:
    Plasma5();
    Plasma5(const int idx, const sol::table &def);
    ~Plasma5();

    enum class TimeIntegrator {
        ForwardsEuler=0,
        BackwardsEuler,
    };

    virtual void solve(MFP* mfp, const Real dt) const override;

    void explicit_solve(MFP* mfp, const Real dt) const;
    void implicit_solve(MFP* mfp, const Real dt) const;

    virtual std::string get_tag() const override {return tag;}

    FieldState* field;
    Vector<HydroState*> species;

    static TimeIntegrator solver;

    static bool registered;
    static std::string tag;
};

#endif // PLASMA5_H
