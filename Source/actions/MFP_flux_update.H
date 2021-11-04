#ifndef FLUX_UPDATE_H
#define FLUX_UPDATE_H

#include "MFP_action.H"
#include "MFP_eulerian.H"

class FluxUpdate : public Action
{
public:
    FluxUpdate();
    FluxUpdate(const int idx, const sol::table &def);
    ~FluxUpdate();

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::Flux;}

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_spatial_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt, const Real flux_register_scale) override;

    Vector<EulerianState*> states;


    static bool registered;
    static std::string tag;
};

#endif // CTU_H
