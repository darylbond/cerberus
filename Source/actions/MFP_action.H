#ifndef ACTION_H
#define ACTION_H

#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include "sol.hpp"
#include "MFP_factory.H"

class Action;
class MFP;

using namespace amrex;

class Action
{
public:
    Action();
    ~Action();

    int src_idx;
    std::string name;
    Vector<int> state_indexes;

    // calculate time and space derivatives for situations where we can accumulate the effects
    virtual void calc_time_derivative(MFP* mfp, Vector<std::pair<int,MultiFab>>& dU, const Real time, const Real dt) const {}
    virtual void calc_spatial_derivative(MFP* mfp, Vector<std::pair<int,MultiFab>>& dU, const Real time, const Real dt) const {}

    // directly update quantities where an accumulation is not possible, these will overwrite stuff!!
    virtual void apply_time_derivative(MFP* mfp, const Real time, const Real dt) const {}
    virtual void apply_spatial_derivative(MFP* mfp, const Real time, const Real dt) const {}



    virtual std::string get_tag() const = 0;
};

template <typename D>
std::unique_ptr<Action> ActionBuilder(const sol::table& def)
{
    int idx = def["src_idx"];

    if (def["type"] == D::tag) {
        return std::unique_ptr<D>(new D(idx, def));
    } else {
        return nullptr;
    }
}

ClassFactory<Action> &GetActionFactory();

#endif // ACTION_H
