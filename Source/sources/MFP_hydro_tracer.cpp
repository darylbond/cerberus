#include "MFP_hydro_tracer.H"
#include "sol.hpp"

std::string HydroTracer::tag = "hydro_tracer";
bool HydroTracer::registered = GetSourceFactory().Register(HydroTracer::tag, SourceBuilder<HydroTracer>);


HydroTracer::HydroTracer(){}
HydroTracer::~HydroTracer(){}

HydroTracer::HydroTracer(const int idx, const sol::table &def)
{
    src_idx = idx;
    name = def["name"];

    std::string hydro_name = def["fluid"];
    std::string particle_name = def["particles"];

    for (const auto& i : MFP::eulerian_states) {
        EulerianState& istate = EulerianState::get_state_global(i);
        if (istate.get_type() == State::StateType::Hydro) {
            if (istate.name == hydro_name) {
                hydro_state = static_cast<HydroState*>(&istate);
            }
        }
    }

    for (const auto& i : MFP::lagrangian_states) {
        LagrangianState& istate = LagrangianState::get_state_global(i);
        if (istate.get_type() == State::StateType::TracerParticle) {
            if (istate.name == particle_name) {
                tracer_state = static_cast<TracerParticle*>(&istate);
            }
        }
    }
}


void HydroTracer::push_particles(TParTileType& ptile,
                                 const FArrayBox& prim,
                                 const Geometry geom,
                                 const Real dt
                                 #ifdef AMREX_USE_EB
                                 ,const EBCellFlagFab& flag
                                 #endif
                                 ) const
{
    BL_PROFILE("MFP::push_particles");

    // advance particle locations

//    const Real          strttime = amrex::second();
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();

    for (int ipass = 0; ipass < 2; ipass++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif


        auto& aos  = ptile.GetArrayOfStructs();
        const int n          = aos.numParticles();


        const auto p4 = prim.array();
        auto  p_pbox = aos().data();

#ifdef AMREX_USE_EB

        const auto f4 = flag.array();
#endif
        amrex::ParallelFor(n,
                           [=] AMREX_GPU_DEVICE (int i)
        {
            TParticleType& p  = p_pbox[i];
            if (p.id() <= 0) return;
            Real v[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0,0.0,0.0)};

            // implement particle boundary conditions here??

            // calculate where we are in index space and cell local space
            Array<Real,3> loc = {0.0,0.0,0.0};
            Array<int,3> iloc = {0,0,0};

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                loc[d] = (p.pos(d) - plo[d]) * dxi[d] - 0.5;
                iloc[d] = static_cast<int>(amrex::Math::floor(loc[d]));
                loc[d] -= iloc[d];
            }
#ifdef AMREX_USE_EB
            if (f4(iloc[0], iloc[1], iloc[2]).isCovered()) {
                p.id() = -p.id();
                return;
            }
#endif

            for (int d1=0; d1<AMREX_SPACEDIM; ++d1) {
                v[d1] = p4(iloc[0], iloc[1], iloc[2],d1); // cell centre value
            }

            // apply updates to particle position and velocity

            if (ipass == 0) {
                for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                    p.rdata(dim) = p.pos(dim);
                    p.pos(dim) += 0.5*dt*v[dim];
                }
            } else  {
                for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                    p.rdata(dim) = p.rdata(dim) + dt*v[dim];
                    p.rdata(dim) = v[dim];
                }
            }
        });
    }

    return;
}






void HydroTracer::solve(MFP* mfp, const Real dt) const
{
    BL_PROFILE("HydroTracer::solve");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    MultiFab& hydro_data = mfp->get_new_data(hydro_state->data_idx);


    FArrayBox vel;

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
        EBData& eb = mfp->get_eb_data(hydro_state->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        const EBCellFlagFab& flag = eb.flags[mfi];
        if (vfrac.getType() == FabType::covered) continue;
#endif

        hydro_state->calc_velocity(box,
                                   hydro_data[mfi],
                                   vel,
                           #ifdef AMREX_USE_EB
                                   vfrac
                           #endif
                                   );

        TParTileType& ptile = tracer_state->particles->ParticlesAt(mfp->get_level(), mfi);

        push_particles(ptile,
                       vel,
                       mfp->Geom(),
                       dt
               #ifdef AMREX_USE_EB
                       ,flag
               #endif
                       );




                // update the cost function
                wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
