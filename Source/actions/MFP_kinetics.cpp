#ifdef EILMER_GAS
#include "MFP_kinetics.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"

std::string GasKinetics::tag = "gas";
bool GasKinetics::registered = GetActionFactory().Register(GasKinetics::tag, ActionBuilder<GasKinetics>);


// D runtime
extern "C" int cwrap_gas_init();
extern "C" int gas_model_new(char* file_name);
extern "C" int gas_model_n_species(int gm_i);
extern "C" int thermochemical_reactor_new(int gm_i, char* filename1, char* filename2);


GasKinetics::GasKinetics(){}
GasKinetics::~GasKinetics(){}

GasKinetics::GasKinetics(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    const sol::table state_names = def["states"];

    neutrals = &HydroState::get_state(state_names["neutral"]);
    if (neutrals->get_type() != State::StateType::Hydro)
        Abort("An invalid neutral state has been defined for the 'gas' action "+name);
    state_indexes.push_back(neutrals->global_idx);
    neutrals->associated_actions.push_back(action_idx);

    ions = &HydroState::get_state(state_names["ion"]);
    if (ions->get_type() != State::StateType::Hydro)
        Abort("An invalid ion state has been defined for the 'gas' action "+name);
    state_indexes.push_back(ions->global_idx);
    ions->associated_actions.push_back(action_idx);

    electrons = &HydroState::get_state(state_names["electron"]);
    if (electrons->get_type() != State::StateType::Hydro)
        Abort("An invalid electron state has been defined for the 'gas' action "+name);
    state_indexes.push_back(electrons->global_idx);
    electrons->associated_actions.push_back(action_idx);


    // gas construction

    cwrap_gas_init();

    std::string gas_model_file = def.get_or<std::string>("gas_model","");
    if (gas_model_file.empty()) Abort("Action '"+name+"' requires 'gas_model' to be defined (a lua file)");

    gas_model_id = gas_model_new(gas_model_file.data());

    if (gas_model_id < 0) Abort("Action '"+name+"' has failed when trying to create a new gas model");

    n_species = gas_model_n_species(gas_model_id);

    if (neutrals->n_species != n_species)
        Abort("Number of species in '"+neutrals->name+"' doesn't match the gas model in '"+name+"'");


    std::string chemistry_update_file = def.get_or<std::string>("chemistry_update","");
    if (chemistry_update_file.empty()) Abort("Action '"+name+"' requires 'chemistry_update' to be defined (a lua file)");

    std::string chemistry_update_file2 = def.get_or<std::string>("chemistry_update2","");

    thermochemical_reactor_id = thermochemical_reactor_new(gas_model_id, chemistry_update_file.data(), chemistry_update_file2.data());

    if (thermochemical_reactor_id < 0) Abort("Action '"+name+"' has failed when trying to create a thermochemical reactor");


    return;
}

void GasKinetics::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("Collisions::get_data");

    Vector<Array<int,2>> options = {
        {neutrals->global_idx, 0},
        {ions->global_idx, 0},
        {electrons->global_idx, 0}
    };

    Action::get_data(mfp, options, update, time);

}


void GasKinetics::calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("Collisions::calc_time_derivative");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);


    update[neutrals->data_idx].dU_status = UpdateData::Status::Changed;
    update[ions->data_idx].dU_status = UpdateData::Status::Changed;
    update[electrons->data_idx].dU_status = UpdateData::Status::Changed;

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(neutrals->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        const Array4<const Real>& n_U4 = update[neutrals->data_idx].U.array(mfi);
        const Array4<Real>& n_dU4 = update[neutrals->data_idx].dU.array(mfi);

        const Array4<const Real>& i_U4 = update[ions->data_idx].U.array(mfi);
        const Array4<Real>& i_dU4 = update[ions->data_idx].dU.array(mfi);

        const Array4<const Real>& e_U4 = update[electrons->data_idx].U.array(mfi);
        const Array4<Real>& e_dU4 = update[electrons->data_idx].dU.array(mfi);


        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif


                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
#endif
