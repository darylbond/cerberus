#ifdef EILMER_GAS
#ifndef KINETICS_H
#define KINETICS_H

#include "MFP_action.H"
#include "MFP_hydro.H"

class GasKinetics : public Action
{
public:
    GasKinetics();
    GasKinetics(const int idx, const sol::table &def);
    ~GasKinetics();

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_time_derivative(MFP* mfp, Vector<UpdateData>& update,  const Real time, const Real dt) override;

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::GasKinetics;}

    HydroState *neutrals, *ions, *electrons;

    int gas_model_id, thermochemical_reactor_id;
    int n_species;

    static bool registered;
    static std::string tag;
};

#endif // KINETICS_H
#endif
