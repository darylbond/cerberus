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

    virtual void solve(MFP* mfp, const Real dt) const override;

    virtual std::string get_tag() const override {return tag;}

    FieldState* field;
    Vector<HydroState*> species;

    static bool registered;
    static std::string tag;
};

#endif // PLASMA5_H
