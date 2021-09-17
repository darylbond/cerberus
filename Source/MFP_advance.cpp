#include "MFP.H"
#include "MFP_eulerian.H"
#include "MFP_diagnostics.H"

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
        EulerianState &istate = EulerianState::get_state(data_idx);

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

            EulerianState &istate = EulerianState::get_state(data_idx);

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
            active[data_idx] = FabType::regular;
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

//            plot_FAB_1d(prim, "prim-"+num2str(level), true);
//                        plot_FAB_2d(prim, 0, "prim 0", false, true);


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
                               *conserved[idx],
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
        if (CTU) {
            for (int idx=0; idx<n_eulerian; ++idx) {


                if (active[idx] == FabType::covered) continue;


                EulerianState &istate = EulerianState::get_state(idx);

                // recalculate the fluxes
                istate.calc_fluxes(box,
                                   *conserved[idx],
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

            istate.calc_viscous_fluxes(grow(box,num_grow_eb),
                                       fluxes[idx],
                                       pbox, primitives,
                           #ifdef AMREX_USE_EB
                                       *fab_flags[idx],
                           #endif
                                       dx);


            // given all of the fluxes calculate the update to the cell centred values
            const int as_crse = (fr_as_crse[idx] != nullptr);
            const int as_fine = (fr_as_fine[idx] != nullptr);

            nu = local_new[idx].nComp();

            FArrayBox& dU = update[idx][mfi];


#ifdef AMREX_USE_EB

            Array<const FArrayBox*,AMREX_SPACEDIM> afrac, fcent;

            if (active[idx] != FabType::regular) {

                const FArrayBox &bcent = (*istate.eb_data.bndrycent)[mfi];
                const FArrayBox &bnorm = (*istate.eb_data.bndrynorm)[mfi];

                CutFab& bc_idx = istate.eb_data.bndryidx[mfi];

                afrac = {AMREX_D_DECL(&(*istate.eb_data.areafrac[0])[mfi], &(*istate.eb_data.areafrac[1])[mfi], &(*istate.eb_data.areafrac[2])[mfi])};


                // calculate the flux through cut cell faces

                istate.calc_wall_fluxes(box,
                                        primitives[idx],
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

                fcent = {AMREX_D_DECL(&(*istate.eb_data.facecent[0])[mfi], &(*istate.eb_data.facecent[1])[mfi], &(*m_istate.eb_data.facecent[2])[mfi])};

                istate.calc_eb_divergence(box,
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

                istate.merge_cells(box,
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

//            plot_FAB_1d(s2, "cons-"+num2str(level), true);

        }

        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }

    return;
}
