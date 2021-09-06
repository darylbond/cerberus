#include "MFP.H"
#include "MFP_state.H"



using namespace amrex;

sol::state MFP::lua;
std::string MFP::lua_script;

Real MFP::force_dt = 0.0;
Real MFP::cfl = 0.0;

int MFP::Cost_Idx;
bool MFP::archive_checkpoint = true;
int MFP::verbosity;
int MFP::linear_solver_verbosity;

Vector<MFP::RefineBoxType> MFP::refine_box_type;
Vector<amrex::RealBox> MFP::refine_boxes;
bool MFP::only_refine_in_box;
bool MFP::derefine_box;

bool MFP::refine_cutcells;

bool MFP::zero_dimensional = false;

IntVect MFP::tile_size;

std::map<std::string,Array<int,AMREX_SPACEDIM+1>> MFP::plot_variables;
Vector<std::pair<std::string,Optional3D1VFunction>> MFP::plot_functions;

constexpr int MFP::level_mask_interior;
constexpr int MFP::level_mask_covered;
constexpr int MFP::level_mask_notcovered;
constexpr int MFP::level_mask_physbnd;

Vector<std::unique_ptr<State>> MFP::states;
Vector<std::string> MFP::state_names;
std::map<std::string, int> MFP::state_index;


#ifdef AMREX_USE_EB
Vector<DataEB> MFP::eb_def;
#endif

Real MFP::x_ref = 1.0;
Real MFP::n_ref = 1.0;
Real MFP::m_ref = 1.0;
Real MFP::rho_ref = 1.0;
Real MFP::T_ref = 1.0;
Real MFP::u_ref = 1.0;
Real MFP::n0 = 1.0;
Real MFP::lightspeed = 1.0;
Real MFP::beta = 1.0;
Real MFP::skin_depth = 1.0;
Real MFP::Larmor = 1.0;
Real MFP::Debye = 1.0;

bool MFP::plasma_params_set = false;

MFP::MFP() {}

MFP::MFP(Amr& papa, int lev, const Geometry& level_geom, const BoxArray& bl,
         const DistributionMapping& dm, Real time)
    : AmrLevel(papa, lev, level_geom, bl, dm, time)
{
    if (level > 0) {

        flux_reg.resize(get_num_fluid_states());

        for (const auto& state : states) {
            int idx = state->data_idx;
            if (idx >= 0) {
                flux_reg[idx].define(
                            bl, papa.boxArray(level - 1),
                            dm, papa.DistributionMap(level - 1),
                            level_geom, papa.Geom(level - 1),
                            papa.refRatio(level - 1), level,
                            desc_lst[idx].nComp());
            }
        }

    }

    build_eb();
}

MFP::~MFP(){}

void MFP::init(AmrLevel& old)
{
    auto& oldlev = dynamic_cast<MFP&>(old);

    Real dt_new = parent->dtLevel(level);
    Real cur_time = oldlev.state[0].curTime();
    Real prev_time = oldlev.state[0].prevTime();
    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    MultiFab& C_new = get_new_data(Cost_Idx);
    FillPatch(old, C_new, 0, cur_time, Cost_Idx, 0, 1);

    for (const auto& state : states) {
        int idx = state->data_idx;
        if (idx >= 0) {
            MultiFab& S_new = get_new_data(idx);
            FillPatch(old, S_new, 0, cur_time, idx, 0, S_new.nComp());
        }
    }
}

void MFP::init()
{
    Real dt = parent->dtLevel(level);
    Real cur_time = getLevel(level - 1).state[0].curTime();
    Real prev_time = getLevel(level - 1).state[0].prevTime();
    Real dt_old = (cur_time - prev_time) /
                  static_cast<Real>(parent->MaxRefRatio(level - 1));
    setTimeLevel(cur_time, dt_old, dt);

    MultiFab& C_new = get_new_data(Cost_Idx);
    FillCoarsePatch(C_new, 0, cur_time, Cost_Idx, 0, 1);

    for (const auto& state : states) {
        int idx = state->data_idx;
        if (idx >= 0) {
            MultiFab& S_new = get_new_data(idx);
            FillCoarsePatch(S_new, 0, cur_time, idx, 0, S_new.nComp());
        }
    }
}

void MFP::initData()
{
    BL_PROFILE("MFP::initData()");

    for (const auto& state : states) {

        state->init_data(this);

    }

    // initialise cost state
    MultiFab& C_new = get_new_data(Cost_Idx);
    C_new.setVal(1.0);

}

Real MFP::advance(Real time, Real dt, int iteration, int ncycle)
{
    return dt;
}

Real MFP::estTimeStep() {
    BL_PROFILE("MFP::estTimeStep()");

    if (force_dt > 0.0) {
        return force_dt;
    }

    Real estdt = std::numeric_limits<Real>::max();

    // handle all of the intrinsic time step limitations of a state
    for (auto& state : states) {
        estdt = std::min(estdt, state->get_allowed_time_step(this));
    }

    // do source term limitations here...

    estdt *= cfl;
    ParallelDescriptor::ReduceRealMin(estdt);

    return estdt;
}

void MFP::computeInitialDt(int finest_level, int sub_cycle,
                           Vector<int>& n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>& dt_level, Real stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0) {
        return;
    }

    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        dt_level[i] = getLevel(i).estTimeStep();
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0, n_factor * dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
    Real cur_time = state[0].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps)) dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void MFP::computeNewDt(int finest_level, int sub_cycle, Vector<int>& n_cycle,
                       const Vector<IntVect>& ref_ratio, Vector<Real>& dt_min,
                       Vector<Real>& dt_level, Real stop_time,
                       int post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0) {
        return;
    }

    for (int i = 0; i <= finest_level; i++) {
        dt_min[i] = getLevel(i).estTimeStep();
    }

    if (post_regrid_flag == 1) {
        //
        // Limit dt's by pre-regrid dt
        //
        for (int i = 0; i <= finest_level; i++) {
            dt_min[i] = std::min(dt_min[i], dt_level[i]);
        }
    } else {
        //
        // Limit dt's by change_max * old dt
        //
        static Real change_max = 1.1;
        for (int i = 0; i <= finest_level; i++) {
            dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
        }
    }

    //
    // Find the minimum over all levels
    //
    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0, n_factor * dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
    Real cur_time = state[0].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void MFP::post_regrid(int lbase, int new_finest)
{


}

void MFP::post_timestep(int iteration)
{
    // reflux
    if (level < parent->finestLevel()) {
        MFP& fine_level = getLevel(level + 1);


        for (const auto& state : states) {
            int idx = state->data_idx;
            State& fine_state = fine_level.get_state(state->global_idx);
            if (idx >= 0) {
                MultiFab& S_crse = get_new_data(idx);
#ifdef AMREX_USE_EB
                MultiFab& S_fine = fine_level.get_new_data(idx);
                fine_level.flux_reg[idx].Reflux(S_crse, state->eb_data.volfrac, S_fine,
                                                fine_state.eb_data.volfrac);
#else
                fine_level.flux_reg[idx].Reflux(S_crse, 0);
#endif


            }
        }
    }

    if (level < parent->finestLevel()) {
        avgDown();
    }

}

void MFP::avgDown()
{
    BL_PROFILE("MFP::avgDown()");

    if (level == parent->finestLevel()) return;

    auto& fine_lev = getLevel(level + 1);

    for (const auto& state : states) {
        int idx = state->data_idx;
        MultiFab& S_crse = get_new_data(idx);
        MultiFab& S_fine = fine_lev.get_new_data(idx);

#ifdef AMREX_USE_EB
        MultiFab volume(S_fine.boxArray(), S_fine.DistributionMap(), 1, 0);
        volume.setVal(1.0);

        State& fine_state = fine_lev.get_state(state->global_idx);

        amrex::EB_average_down(S_fine, S_crse, volume, fine_state.eb_data.volfrac, 0,
                               S_fine.nComp(), fine_ratio);
#else
        amrex::average_down(S_fine, S_crse, 0, S_fine.nComp(), fine_ratio);
#endif

    }
}

void MFP::postCoarseTimeStep(Real time)
{


}

void MFP::post_init(Real)
{


}

void MFP::post_restart()
{

}

void MFP::errorEst(TagBoxArray& tags, int, int, Real time, int, int)
{
    BL_PROFILE("MFP::errorEst()");
}


State& MFP::get_state(const std::string& name) {
    if ( state_index.find(name) == state_index.end() ) {
        Abort("Attempting to reference a state that doesn't exist");
    }
    return *states[state_index[name]];
}

int MFP::get_num_fluid_states()
{
    int n = 0;

    for (const auto& state: states) {
        if (state->data_idx > 0) {
            n += 1;
        }
    }

    return n;
}
