#ifdef EILMER_GAS
#include "MFP_kinetics.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"
#include <algorithm>

std::string GasKinetics::tag = "gas";
bool GasKinetics::registered = GetActionFactory().Register(GasKinetics::tag, ActionBuilder<GasKinetics>);


// D runtime
extern "C" int cwrap_gas_init();
extern "C" int gas_model_new(const char* file_name);
extern "C" int gas_model_n_species(int gm_i);
extern "C" int gas_model_species_name(int gm_i, int isp, char* dest_name, int* n);
extern "C" int gas_state_new(int gm_i);
extern "C" int thermochemical_reactor_new(int gm_i, const char* filename1, const char* filename2);
extern "C" int gas_state_set_scalar_field(int gs_i, const char* field_name, double value);
extern "C" int gas_state_get_scalar_field(int gs_i, const char* field_name, double* value);
extern "C" int gas_state_set_array_field(int gs_i, const char* field_name, double* values, int n);
extern "C" int gas_state_get_array_field(int gs_i, const char* field_name, double* values, int n);
extern "C" int gas_model_gas_state_update_thermo_from_rhou(int gm_i, int gs_i);
extern "C" int thermochemical_reactor_gas_state_update(int cr_i, int gs_i, double t_interval, double* dt_suggest);

GasKinetics::GasKinetics(){}
GasKinetics::~GasKinetics(){}

GasKinetics::GasKinetics(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    const sol::table state_names = def["states"];

    for (const auto& key_value_pair : state_names) {
        std::string state_name = key_value_pair.second.as<std::string>();
        State& istate = MFP::get_state(state_name);

        switch (istate.get_type()) {
        case State::StateType::Hydro:
            states.push_back(static_cast<HydroState*>(&istate));
            state_indexes.push_back(istate.global_idx);
            istate.associated_actions.push_back(action_idx);
            break;
        default:
            Abort("An invalid state has been defined for the 'gas' action '"+name+"'");
        }
    }

    num_states = states.size();

    // gas construction

    cwrap_gas_init();

    std::string gas_model_file = def.get_or<std::string>("gas_model","");
    if (gas_model_file.empty()) Abort("Action '"+name+"' requires 'gas_model' to be defined (a lua file)");

    gas_model_id = gas_model_new(gas_model_file.data());

    if (gas_model_id < 0) Abort("Action '"+name+"' has failed when trying to create a new gas model");

    n_species = gas_model_n_species(gas_model_id);

    species_info.resize(n_species);

    std::pair<bool, int > found;
    int n;
    for (int i=0; i<n_species; ++i) {
        std::string& sname = species_info[i].name;
        sname.resize(10);
        gas_model_species_name(gas_model_id, i, sname.data(), &n);
        sname.resize(n);

        // search through the neutrals and ions to get where it is located (neutral/ion, idx)
        for (int s=0; s<num_states; ++s) {
            found = findInVector(states[s]->comp_names, sname);
            if (found.first) {
                species_info[i].species_idx = s;
                species_info[i].alpha_idx = found.second;
                break;
            }
        }

        if (!found.first)
            Abort("Couldn't find sub-component for '"+sname+"'");
    }

    gas_state_id = gas_state_new(gas_model_id);

    if (gas_state_id < 0) Abort("Action '"+name+"' has failed when trying to create a new gas state");


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

    Vector<Array<int,2>> options(states.size());

    for (size_t s=0; s<num_states; ++s) {
        options.push_back({states[s]->global_idx, 0});
    }

    Action::get_data(mfp, options, update, time);

}


void GasKinetics::calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("Collisions::calc_time_derivative");

    const Real* dx = mfp->Geom().CellSize();

    Real volume = AMREX_D_TERM(dx[0]*MFP::x_ref,*dx[1]*MFP::x_ref,*dx[2]*MFP::x_ref);

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    Vector<Array4<const Real>> U4(num_states);
    Vector<Array4<Real>> dU4(num_states);
    Vector<Vector<Real>> U(num_states), alpha(num_states), accounting(num_states);

    for (int si=0; si<num_states; ++si) {
        update[states[si]->data_idx].dU_status = UpdateData::Status::Changed;
        U[si].resize(states[si]->n_cons());
        alpha[si].resize(states[si]->n_species);
        accounting[si].resize(states[si]->n_species);
    }

    // overall gas state
    Real rho_sum, nrg_sum; // note that u is the specific internal energy
    Vector<Real> massf(n_species);

    Vector<Real> s_massf(num_states), s_T(num_states);

    Real dt_suggest;

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(states[0]->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif
        for (int si=0; si<num_states; ++si) {
            U4[si] = update[states[si]->data_idx].U.array(mfi);
            dU4[si] = update[states[si]->data_idx].dU.array(mfi);
        }

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) continue;
#endif

                    // get the conserved quantities and the alpha values
                    for (int si=0; si<num_states; ++si) {
                        Vector<Real>& UU = U[si];
                        const Array4<const Real>& UU4 = U4[si];
                        for (size_t l=0; l<UU.size(); ++l) {
                            UU[l] = UU4(i,j,k,l);
                        }

                        if (UU[+HydroDef::ConsIdx::Density] < states[si]->effective_zero) {
                            std::fill(alpha[si].begin(), alpha[si].end(), 0.0);
                            s_T[si] = 0.0;
                        } else {
                            states[si]->get_alpha_fractions_from_cons(UU, alpha[si]);
                            s_T[si] = states[si]->get_temperature_from_cons(UU);
                        }
                    }

                    // get the overall sum of mass and energy held by the gas state
                    nrg_sum = 0.0;
                    std::fill(s_massf.begin(), s_massf.end(), 0.0);
                    for (int isp=0; isp<n_species; ++isp) {
                        SpeciesInfo& si = species_info[isp];
                        s_massf[si.species_idx] += alpha[si.species_idx][si.alpha_idx]*U[si.species_idx][+HydroDef::ConsIdx::Density];
                        nrg_sum                 += alpha[si.species_idx][si.alpha_idx]*U[si.species_idx][+HydroDef::ConsIdx::Eden];
                    }

                    rho_sum = sum_vec(s_massf);

                    // get the mass fractions
                    for (int isp=0; isp<n_species; ++isp) {
                        SpeciesInfo& si = species_info[isp];
                        massf[isp] = alpha[si.species_idx][si.alpha_idx]*U[si.species_idx][+HydroDef::ConsIdx::Density]/rho_sum;
                    }

                    // calculate the mass weighted temperature
                    for (int si=0; si<num_states; ++si) {
                        s_massf[si] /= rho_sum;
                        s_T[si] *= s_massf[si];
                    }
                    const Real guess_T = sum_vec(s_T);



                    // set properties - note the conversion to dimensional values
                    gas_state_set_scalar_field(gas_model_id, "rho", rho_sum*MFP::rho_ref);
                    gas_state_set_scalar_field(gas_model_id, "u", nrg_sum*MFP::prs_ref/(rho_sum*MFP::rho_ref));
                    gas_state_set_scalar_field(gas_model_id, "T", guess_T*MFP::T_ref); // need an initial guess at the temperature
                    gas_state_set_array_field(gas_state_id, "massf", massf.data(), n_species);

                    // update model and solve for update
                    gas_model_gas_state_update_thermo_from_rhou(gas_model_id, gas_state_id);
                    thermochemical_reactor_gas_state_update(thermochemical_reactor_id, gas_state_id, dt*MFP::t_ref, &dt_suggest);

                    // retrieve mass fraction data
                    gas_state_get_array_field(gas_state_id, "massf", massf.data(), n_species);

                    // accounting for changes in mass and tracers (mass fractions)
                    for (int si=0; si<num_states; ++si) {
                        for (int ci=0; ci<alpha[si].size(); ++ci) {
                            accounting[si][ci] = U[si][+HydroDef::ConsIdx::Density]*alpha[si][ci]; // original mass
                        }
                    }

                    for (int isp=0; isp<n_species; ++isp) {
                        SpeciesInfo& si = species_info[isp];
                        accounting[si.species_idx][si.alpha_idx] = massf[isp]*rho_sum;
                    }

                    for (int si=0; si<num_states; ++si) {
                        const Real rho = sum_vec(accounting[si]);
                        dU4[si](i,j,k,+HydroDef::ConsIdx::Density) += rho - U4[si](i,j,k,+HydroDef::ConsIdx::Density);
                        for (int ci=0; ci<states[si]->n_tracers; ++ci) {
                            dU4[si](i,j,k,+HydroDef::ConsIdx::NUM + ci) += accounting[si][ci]  - U4[si](i,j,k,+HydroDef::ConsIdx::NUM + ci);
                        }
                    }

                    // accounting for changes in energy density
                    for (int si=0; si<num_states; ++si) {
                        for (int ci=0; ci<alpha[si].size(); ++ci) {
                            accounting[si][ci] = U[si][+HydroDef::ConsIdx::Eden]*alpha[si][ci]; // original energy density
                        }
                    }

                    for (int isp=0; isp<n_species; ++isp) {
                        SpeciesInfo& si = species_info[isp];
                        accounting[si.species_idx][si.alpha_idx] = massf[isp]*nrg_sum;
                    }

                    for (int si=0; si<num_states; ++si) {
                        const Real u = sum_vec(accounting[si]);
                        dU4[si](i,j,k,+HydroDef::ConsIdx::Eden) += u - U4[si](i,j,k,+HydroDef::ConsIdx::Eden);
                    }


                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
#endif
