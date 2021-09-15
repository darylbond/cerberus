#include "MFP.H"
#include "MFP_state.H"
#include "MFP_eulerian.H"


using namespace amrex;

MFP::TimeIntegrator MFP::time_integration_scheme;

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
Vector<size_t> MFP::eulerian_states;
Vector<size_t> MFP::lagrangian_states;


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
    BL_PROFILE("MFP::advance()");

    MultiFab& C_new = get_new_data(Cost_Idx);
    C_new.setVal(0.0);

    switch (time_integration_scheme) {
    case TimeIntegrator::Euler:
        advance_one_step(time, dt, iteration, ncycle);
        break;
    case TimeIntegrator::CTU:
        advance_one_step(time, dt, iteration, ncycle, true);
        break;
    default:
        advance_one_step(time, dt, iteration, ncycle);
    }

    for (int i = 0; i < states.size(); ++i) {
        state[i].setNewTimeLevel(time+dt);
    }


    // do some lua garbage collection
    //    lua.script("collectgarbage('collect')");
    //    lua.script("collectgarbage('collect')");

    return dt;
}

void MFP::advance_one_step(Real time, Real dt, int iteration, int ncycle, bool CTU)
{

    BL_PROFILE("MFP::advance_euler");

#ifdef AMREX_USE_EB
    constexpr int num_grow_eb = 2;
#else
    constexpr int num_grow_eb = 0;
#endif

    const size_t n_states = states.size();
    const size_t n_eulerian = eulerian_states.size();

    // get the maximum wave speed from the time step and cell spacing
    RealVect speed_max;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        speed_max[d] = dt/(cfl*geom.CellSize(d));
    }

    // ==========================================================================
    // all of level storage

    Vector<MultiFab*> global_new(n_eulerian, nullptr);
    Vector<MultiFab> local_new(n_eulerian); // local 'new', gathered from whatever is most recent from the state vector
    Vector<MultiFab> update(n_eulerian); // fluxes, sources
    Vector<FabType> active(n_eulerian); // flag for calculation

#ifdef AMREX_USE_EB
    Vector<EBFluxRegister*> fr_as_crse(n_eulerian, nullptr);
    Vector<EBFluxRegister*> fr_as_fine(n_eulerian, nullptr);
#else
    Vector<YAFluxRegister*> fr_as_crse(n_eulerian, nullptr);
    Vector<YAFluxRegister*> fr_as_fine(n_eulerian, nullptr);
#endif

    for (int data_idx = 0; data_idx < n_eulerian; ++data_idx) {
        int global_idx = eulerian_states[data_idx];

        EulerianState &istate = EulerianState::get_state(global_idx);

        int ns = desc_lst[data_idx].nComp();
        int ng = istate.get_num_grow() + num_grow_eb;

        local_new[data_idx].define(grids, dmap, ns, ng, MFInfo(),Factory());
        update[data_idx].define(grids, dmap, ns, num_grow_eb, MFInfo(),Factory());
        global_new[data_idx] = &get_new_data(data_idx);

        // get a full array of data at this level
#ifdef AMREX_USE_EB
        EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate.eb2_index));
#endif
        FillPatch(*this, local_new[data_idx], ng, time, data_idx, 0, ns);

        //        plot_FAB_2d(local_new[idx], 0, 0, "cons density", false, false);
        //        plot_FAB_2d(flag, "flag", true);

        if (istate.reflux && level < parent->finestLevel()) {
            MFP& fine_level = getLevel(level + 1);
            fr_as_crse[data_idx] = &fine_level.flux_reg[data_idx];
        }

        if (istate.reflux && level > 0) {
            fr_as_fine[data_idx] = &flux_reg[data_idx];
        }

        if (fr_as_crse[data_idx]) {
            fr_as_crse[data_idx]->reset();
        }
    }

    // ==========================================================================
    // per state storage

    Vector<FArrayBox*> conserved(n_eulerian);
    Vector<FArrayBox> primitives(n_eulerian);
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> R_lo(n_eulerian), R_hi(n_eulerian);
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> fluxes(n_eulerian);

    int nu; // per state size

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& cost = get_new_data(Cost_Idx);

#ifdef AMREX_USE_EB
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> wall_fluxes(n_eulerian);
    Vector<const EBCellFlagFab*> fab_flags(n_eulerian);
    Vector<const FArrayBox*> fab_vfrac(n_eulerian);
#endif

    // iterate over all of the FABs within the level performing reconstruction, flux calculation,
    // and updating the cell-centred data according to the resulting divergence

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();

        // region over which to perform reconstruction for face values
        const Box rbox = amrex::grow(box, 1 + num_grow_eb);

        // ==========================================================================
        // 1. iterate over all states to set-up the data required for flux calculation
        for (int data_idx=0; data_idx<n_eulerian; ++data_idx) {

            int global_idx = eulerian_states[data_idx];

            EulerianState &istate = EulerianState::get_state(global_idx);

            // get a pointer to the conserved quantities
            conserved[data_idx] = &local_new[data_idx][mfi];

#ifdef AMREX_USE_EB
            // get the EB data required for later calls
            const EBCellFlagFab& flag = istate.eb_data.flags[mfi]; fab_flags[data_idx] = &flag;
            const FArrayBox& vfrac = istate.eb_data.volfrac[mfi]; fab_vfrac[data_idx] = &vfrac;

            //            plot_FAB_2d(flag, "flag", false);
            //            plot_FAB_2d(vfrac, 0, "vfrac", false, true);

            // check if this box is active
            active[data_idx] = flag.getType(rbox);
            if (active[data_idx] == FabType::covered)
                continue;
#else
            active[idx] = FabType::regular;
#endif
            // final check to see if this state actually has any transport
            // this needs to be here (rather than higher up) as we need
            // certain variables to be defined for other uses (face sources)
            if (!istate.is_transported()) {
                active[data_idx] = FabType::covered;
                continue;
            }

            // region over which to get cell centered primitives for reconstruction
            const Box pbox = amrex::grow(box, istate.num_grow + num_grow_eb);

            // ===============================================
            // 1.1 Calculate primitive values within each cell

            FArrayBox& cons = local_new[data_idx][mfi];
            FArrayBox& prim = primitives[data_idx];

            istate.calc_primitives(pbox,
                                   cons,
                                   prim,
                                   dx,
                                   time,
                                   prob_lo
                       #ifdef AMREX_USE_EB
                                   ,vfrac
                       #endif
                                   );

            //            plot_FAB_2d(prim, 0, "prim 0", false, true);


            // fill in any cells that need special boundary values
            istate.update_boundary_cells(pbox,
                                         geom,
                                         prim,
                             #ifdef AMREX_USE_EB
                                         vfrac,
                             #endif
                                         time);



            // =======================================
            // 1.2 Calculate reconstructed face values

            // each cell has a hi and lo side in each direction

            // calculate the reconstructed face values
            istate.calc_reconstruction(rbox,
                                       prim,
                                       R_lo[data_idx],
                                       R_hi[data_idx]
                           #ifdef AMREX_USE_EB
                                       ,flag
                                       ,vfrac
                           #endif
                                       );

            // ===================================================================
            // update the face values to time t+1/2 based on the local wave speeds

            // TODO: this step currently assumes a single speed for all components and
            // should be updated to calculate the correct characteristic speeds
            // for each component individually
            if (CTU) {
                istate.calc_time_averaged_faces(rbox,
                                                prim,
                                                R_lo[data_idx],
                                                R_hi[data_idx],
                                #ifdef AMREX_USE_EB
                                                flag,
                                #endif
                                                dx,
                                                dt);
            }

            // =======================================
            // 1.3 Particle update

            // update particle locations and velocities using the reconstructed
            // face values to interpolate the local velocity for each particle
#ifdef AMREX_PARTICLES
            //            if (gd.do_tracer_particles && istate.particle_index > -1) {
            //                AmrTracerParticleContainer* pc = particles[istate.particle_index];
            //                if (pc) {

            //                    // grab the starting index for velocity
            //                    const int vel_idx = istate.get_prim_vector_idx()[0];

            //                    // grab the tile of particles
            //                    auto& ptile = pc->ParticlesAt(level, mfi);

            //                    // update the position and velocity of the tracer particles
            //                    // using the reconstructed primitive values on cell faces
            //                    push_particles(ptile,
            //                                   prim,
            //                                   R_lo[data_idx],
            //                                   R_hi[data_idx],
            //                                   vel_idx,
            //                                   dt
            //                                   EB_OPTIONAL(,flag)
            //                                   );

            //                }
            //            }
#endif

        }

        // ==========================================================================
        // 2. Apply face value modifications (assumes t = t + dt/2)

        //        if (do_face_src && CTU) {
        //            calc_face_source(rbox,
        //                             conserved,
        //                             R_lo,
        //                             R_hi,
        //                             EB_OPTIONAL(fab_flags,)
        //                             dx,
        //                             time+dt/2,
        //                             dt);
        //        }

        // ==========================================================================
        // 3.1 Setup for flux calculation

        // resize the flux arrays before any get used
        for (int idx=0; idx<n_eulerian; ++idx) {


            if (active[idx] == FabType::covered)
                continue;


            EulerianState& istate = EulerianState::get_state(idx);

            nu = local_new[idx].nComp();
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // get a node centered box in direction d that encloses the original box
                Box fbox = surroundingNodes(box, d);

                if (CTU || (istate.is_viscous() && num_grow_eb > 0)) {
                    // enlarge the fluxes for face corrections
                    // grow in all directions but that of the flux
                    IntVect expand(1);
                    expand.setVal(d, 0);
                    fbox = amrex::grow(fbox,expand);
                }

#ifdef AMREX_USE_EB
                fbox = grow(fbox, num_grow_eb);
                wall_fluxes[idx][d].resize(fbox, nu);
#endif
                fluxes[idx][d].resize(fbox, nu);
            }
        }

        // ==========================================================================
        // 3.2 Calculate fluxes

        //---
        for (int idx=0; idx<n_eulerian; ++idx) {

            if (active[idx] == FabType::covered) {
                update[idx][mfi].setVal(0.0);
                continue;
            }

            EulerianState &istate = EulerianState::get_state(idx);

            istate.update_face_prim(box,
                                    geom,
                                    R_lo[idx],
                                    R_hi[idx],
                        #ifdef AMREX_USE_EB
                                    *fab_flags[idx],
                        #endif
                                    time);
        }

        //---
        for (int idx=0; idx<n_eulerian; ++idx) {


            if (active[idx] == FabType::covered) continue;


            EulerianState &istate = EulerianState::get_state(idx);

            istate.calc_fluxes(box,
                               R_lo[idx],
                               R_hi[idx],
                               fluxes[idx],
                   #ifdef AMREX_USE_EB
                               *fab_flags[idx],
                   #endif
                               dx, dt);

#if AMREX_SPACEDIM > 1
            if (CTU) {
                // now that we have done the fluxes, do we need to correct them???
                // correction is done according to the following steps:
                // 1. calculate the fluxes (above)
                // 2. modify the reconstructed face states using the flux information
                // 3. recalculate the fluxes with the updated face values

                istate.correct_face_prim(grow(box,num_grow_eb),
                                         R_lo[idx],
                                         R_hi[idx],
                                         fluxes[idx],
                         #ifdef AMREX_USE_EB
                                         *fab_flags[idx],
                         #endif
                                         dx, dt);

                // following the update of the face values we need to update any boundary conditions

                istate.update_face_prim(box,
                                        geom,
                                        R_lo[idx],
                                        R_hi[idx],
                        #ifdef AMREX_USE_EB
                                        *fab_flags[idx],
                        #endif
                                        time,
                                        true);
            }
#endif
        }

        //---
#if AMREX_SPACEDIM > 1
        if (do_CTU) {
            for (int idx=0; idx<n_eulerian; ++idx) {


                if (active[idx] == FabType::covered) continue;


                EulerianState &istate = EulerianState::get_state(idx);

                // recalculate the fluxes
                istate.calc_fluxes(box,
                                   R_lo[idx],
                                   R_hi[idx],
                                   fluxes[idx],
                   #ifdef AMREX_USE_EB
                                   *fab_flags[idx],
                   #endif
                                   dx, dt);
            }
        }
#endif

        //---
        for (int idx=0; idx<n_eulerian; ++idx) {


            if (active[idx] == FabType::covered) continue;


            EulerianState &istate = EulerianState::get_state(idx);

            // now calculate any viscous fluxes

            const Box pbox = amrex::grow(box, istate.num_grow + num_grow_eb);
            FArrayBox& prim = primitives[idx];

            istate.calc_viscous_fluxes(grow(box,num_grow_eb),
                                       fluxes[idx],
                                       pbox, primitives,
                           #ifdef AMREX_USE_EB
                                       *fab_flags[idx],
                           #endif
                                       dx);

            // shock tracking
            if (gd.Shock_Idx > 0 && istate.shock_idx > -1) {
                (*mf_shock)[mfi].copy(shock[idx], 0, istate.shock_idx);
            }

            // given all of the fluxes calculate the update to the cell centred values
            const int as_crse = (fr_as_crse[idx] != nullptr);
            const int as_fine = (fr_as_fine[idx] != nullptr);

            nu = local_new[idx].nComp();

            FArrayBox& dU = update[idx][mfi];


#ifdef AMREX_USE_EB

            Array<const FArrayBox*,AMREX_SPACEDIM> afrac, fcent;

            if (active[idx] != FabType::regular) {

                const FArrayBox &bcent = (*getEBData(idx).bndrycent)[mfi];
                const FArrayBox &bnorm = (*getEBData(idx).bndrynorm)[mfi];

                CutFab& bc_idx = getEBData(idx).bndryidx[mfi];

                afrac = {AMREX_D_DECL(&(*getEBData(idx).areafrac[0])[mfi], &(*getEBData(idx).areafrac[1])[mfi], &(*m_getEBData(idx).areafrac[2])[mfi])};


                // calculate the flux through cut cell faces

                istate.calc_wall_fluxes(box,
                                        primitives,
                                        wall_fluxes[idx],
                                        *fab_flags[idx],
                                        bc_idx,
                                        bcent,
                                        bnorm,
                                        afrac,
                                        dx,
                                        dt);

                // calculate divergence, including cut-cells

                FArrayBox dm_as_fine;
                if (as_fine) {
                    dm_as_fine.resize(amrex::grow(box,2),nu);
                    dm_as_fine.setVal(0.0);
                } else {
                    dm_as_fine.resize(Box::TheUnitBox(),nu);
                }

                FArrayBox fab_drho_as_crse(Box::TheUnitBox(),nu);
                FArrayBox* p_drho_as_crse = (fr_as_crse[idx]) ? fr_as_crse[idx]->getCrseData(mfi) : &fab_drho_as_crse;

                IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());
                const IArrayBox* p_rrflag_as_crse = (fr_as_crse[idx]) ? fr_as_crse[idx]->getCrseFlag(mfi) : &fab_rrflag_as_crse;

                fcent = {AMREX_D_DECL(&(*getEBData(idx).facecent[0])[mfi], &(*getEBData(idx).facecent[1])[mfi], &(*m_getEBData(idx).facecent[2])[mfi])};

                istate.eb_div->calc_eb_divergence(box,
                                                  *conserved[idx],
                                                  fluxes[idx],
                                                  wall_fluxes[idx],
                                                  dU,
                                                  *fab_flags[idx],
                                                  *fab_vfrac[idx],
                                                  afrac,
                                                  fcent,
                                                  as_crse,
                                                  as_fine,
                                                  p_rrflag_as_crse,
                                                  level_mask[mfi],
                                                  p_drho_as_crse,
                                                  dm_as_fine,
                                                  dx,
                                                  dt);

                istate.eb_div->merge_cells(box,
                                           *conserved[idx],
                                           dU,
                                           *fab_flags[idx],
                                           *fab_vfrac[idx],
                                           afrac,
                                           as_fine,
                                           dm_as_fine,
                                           level_mask[mfi]);


                if (as_crse) {
                    fr_as_crse[idx]->CrseAdd(mfi,
                    {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])},
                                             dx, dt,
                                             *fab_vfrac[idx],
                                             afrac, RunOn::Cpu);
                }

                if (as_fine) {
                    fr_as_fine[idx]->FineAdd(mfi,
                    {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])},
                                             dx, dt,
                                             *fab_vfrac[idx],
                                             afrac,
                                             dm_as_fine, RunOn::Cpu);
                }

            } else {

#endif
                // calculate divergence

                istate.calc_divergence(box,
                                       fluxes[idx],
                                       dU,
                                       dx,
                                       dt);

                if (fr_as_crse[idx]) {
                    fr_as_crse[idx]->CrseAdd(mfi, {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])}, dx, dt, RunOn::Cpu);
                }

                if (fr_as_fine[idx]) {
                    fr_as_fine[idx]->FineAdd(mfi, {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])}, dx, dt, RunOn::Cpu);
                }

#ifdef AMREX_USE_EB
            }
#endif

        }

        // update the 'new' state held by AmrLevel for all of the states
        // this is performed according to new = new + dt*update

        for (int idx=0; idx<n_eulerian; ++idx) {

            if (active[idx] == FabType::covered) continue;

            FArrayBox& s1 = local_new[idx][mfi];
            FArrayBox& s2 = (*global_new[idx])[mfi];


            s2.linComb(s1, box, 0, update[idx][mfi], box, 0, 1.0, 1.0, box, 0, s1.nComp());

        }

        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }

    return;
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
