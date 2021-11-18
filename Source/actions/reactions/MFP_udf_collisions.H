#ifndef UDF_COLLISIONS_H
#define UDF_COLLISIONS_H

#include "MFP_collisions.H"

class CollisionsUDF : public Collisions
{
public:
    CollisionsUDF();
    CollisionsUDF(const int idx, const sol::table &def);
    ~CollisionsUDF();

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::Collisions;}

    virtual void calc_update() override;

    Vector<Vector<Real> > ccs;

    static bool registered;
    static std::string tag;

    Real n0, Debye, c1, c2, c3, lnC;
};

#endif // UDF_COLLISIONS_H
