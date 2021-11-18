#ifndef ARGON_H
#define ARGON_H

#include "MFP_collisions.H"

class CollisionsArgon : public Collisions
{
public:
    CollisionsArgon();
    CollisionsArgon(const int idx, const sol::table &def);
    ~CollisionsArgon();

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::Collisions;}

    virtual void calc_update();

    Vector<Vector<Real> > ccs;

    static bool registered;
    static std::string tag;

    // [Ar, e-, Ar+]

    enum class SpeciesIndex {
        Ar=0,
        e_minus,
        Ar_plus,
        NUM
    };

    Array<int,+SpeciesIndex::NUM> index;
    Array<Real,+SpeciesIndex::NUM> C, eta, theta, CK, zeta, phi, M, m, d, Z;
};

#endif // UDF_COLLISIONS_H
