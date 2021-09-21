#ifdef AMREX_PARTICLES
#include "MFP_tracer.H"
#include "sol.hpp"

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

Vector<std::string> TracerParticle::particle_real_names = {"x_vel", "y_vel", "z_vel"};
Vector<std::string> TracerParticle::particle_int_names = {};

std::string TracerParticle::tag = "tracer";
bool TracerParticle::registered = GetStateFactory().Register(TracerParticle::tag, StateBuilder<TracerParticle>);

TracerParticle::TracerParticle()
{

}

TracerParticle::TracerParticle(const sol::table& def)
{
    name = def["name"];
    global_idx = def["global_idx"];
    verbosity = def.get_or("verbosity",0);

    // grab the initial positions
    sol::table particle_def = def["particles"];

    for (auto& pd : particle_def) {
        sol::table dat = pd.second;
        initial_positions.push_back({AMREX_D_DECL(dat[1],dat[2],dat[3])});
    }
}

void TracerParticle::init_data(MFP* mfp)
{
    if (mfp->get_level() == 0)
        init(mfp->get_parent());
}

// this should only be called at level 0
void TracerParticle::init(AmrCore* amr_core, bool make_particles)
{
    // generate the particle container
    particles = std::unique_ptr<AmrTParContType>(new AmrTParContType(amr_core));
    particles->SetVerbose(verbosity);

    if (make_particles) {
        // now make the particles

        constexpr int level = 0;

        const Geometry& geom = particles->Geom(level);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        IntVect pos_idx;

        Vector<int> done(initial_positions.size(), 0);

        // iterate over all of the boxes on this level and make particles if they fit into one
        for(MFIter mfi = particles->MakeMFIter(level); mfi.isValid(); ++mfi) {
            Box box = mfi.validbox();

            // get the tile of particles for the local box
            TParTileType& pc = particles->DefineAndReturnParticleTile(level, mfi.index(), mfi.LocalTileIndex());

            for (int pi=0; pi < initial_positions.size(); ++pi) {
                if (done[pi]) continue;

                RealArray& pos = initial_positions[pi];

                // convert position to index
                for (int dim=0; dim<AMREX_SPACEDIM; ++dim) {
                    pos_idx[dim] = std::floor((pos[dim] - plo[dim])*dxi[dim]);
                }

                if (box.contains(pos_idx)) {
                    TParticleType p;
                    p.id()   = TParticleType::NextID();
                    p.cpu()  = ParallelDescriptor::MyProc();
                    AMREX_D_TERM(
                                p.pos(0) = pos[0];,
                            p.pos(1) = pos[1];,
                    p.pos(2) = pos[2];
                    )

                    AMREX_D_TERM(
                                p.rdata(+ParticleIdxR::VX) = 0.0;,
                            p.rdata(+ParticleIdxR::VY) = 0.0;,
                    p.rdata(+ParticleIdxR::VZ) = 0.0;
                    )

                    pc.push_back(p);

                    done[pi] = 1;
                }
            }
        }

        particles->Redistribute();
    }
}

void TracerParticle::checkpoint(const std::string& dir)
{
    particles->Checkpoint(dir, "Particles_"+name, true, particle_real_names, particle_int_names);
}

void TracerParticle::restart(const std::string& dir)
{
    particles->Restart(dir, "Particles_"+name);
}

void TracerParticle::redistribute(int level, int finest_level, int ngrow)
{
    particles->Redistribute(level, finest_level, ngrow);
}

void TracerParticle::clear()
{
    particles.reset();
}


void TracerParticle::push_particles(MFIter& mfi,
                                    const FArrayBox& prim,
                                    Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                    Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                    const int vel_idx,
                                    const Real dt,
                                    const Geometry geom,
                                    const int level
                                    #ifdef AMREX_USE_EB
                                    ,const EBCellFlagFab& flag
                                    #endif
                                    )
{
    BL_PROFILE("TracerParticle::push_particles");

    // advance particle locations

    // grab the tile of particles

    TParTileType& ptile = particles->ParticlesAt(level, mfi);

//    const Real          strttime = amrex::second();
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();

    for (int ipass = 0; ipass < 2; ipass++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif


        auto& aos  = ptile.GetArrayOfStructs();
        const int n  = aos.numParticles();

        // get the left and right cell values
        Array<Array4<const Real>, AMREX_SPACEDIM> rlo4;
        Array<Array4<const Real>, AMREX_SPACEDIM> rhi4;

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            rlo4[d] = rlo[d].array();
            rhi4[d] = rhi[d].array();
        }

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

//            int valid = 0;

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

            // get the velocity at the particle according to the local slopes
            // obtained via the reconstructed cell face values

            for (int d1=0; d1<AMREX_SPACEDIM; ++d1) {
                v[d1] = p4(iloc[0], iloc[1], iloc[2],vel_idx + d1); // cell centre value
                for (int d2=0; d2<AMREX_SPACEDIM; ++d2) {
                    // adjustment due to slopes
                    const Real lo_val = rlo4[d2](iloc[0], iloc[1], iloc[2],vel_idx + d1);
                    const Real hi_val = rhi4[d2](iloc[0], iloc[1], iloc[2],vel_idx + d1);
                    v[d1] += loc[d2]*(hi_val - lo_val);
                }
            }

            // apply updates to particle position and velocity

            if (ipass == 0) {
                for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                    p.rdata(dim) = p.pos(dim);
                    p.pos(dim) += 0.5*dt*v[dim];
                }
            } else  {
                for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                    p.pos(dim) = p.rdata(dim) + dt*v[dim];
                    p.rdata(dim) = v[dim];
                }
            }
        });
    }

    return;
}

void TracerParticle::write_info(nlohmann::json& js) const
{
    LagrangianState::write_info(js);

    js["state_idx"] = state_idx;

}
#endif
