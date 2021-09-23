#ifdef AMREX_PARTICLES
#include "MFP_chargedparticle.H"
#include "MFP_state.H"

Vector<std::string> ChargedParticle::particle_real_names = {"mass", "charge", "x_vel", "y_vel", "z_vel"};
Vector<std::string> ChargedParticle::particle_int_names = {};

std::string ChargedParticle::tag = "charged_particle";
bool ChargedParticle::registered = GetStateFactory().Register(ChargedParticle::tag, StateBuilder<ChargedParticle>);


ChargedParticle::ChargedParticle()
{

}

ChargedParticle::ChargedParticle(const sol::table& def)
{
    name = def["name"];
    global_idx = def["global_idx"];
    verbosity = def.get_or("verbosity",0);

    // grab the initial positions
    sol::table particle_def = def["particles"];

    for (auto& pd : particle_def) {
        sol::table dat = pd.second;
        initial_positions.push_back({AMREX_D_DECL(dat[1],dat[2],dat[3])});
        initial_velocities.push_back({dat[4],dat[5],dat[6]});
    }

    charge = def["charge"];
    mass = def["mass"];
}

void ChargedParticle::init_data(MFP* mfp, const Real time)
{
    if (mfp->get_level() == 0)
        init(mfp->get_parent());
}

// this should only be called at level 0
void ChargedParticle::init(AmrCore* amr_core, bool make_particles)
{
    // generate the particle container
    particles = std::unique_ptr<AmrCParContType>(new AmrCParContType(amr_core));
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
            CParTileType& pc = particles->DefineAndReturnParticleTile(level, mfi.index(), mfi.LocalTileIndex());

            for (int pi=0; pi < initial_positions.size(); ++pi) {
                if (done[pi]) continue;

                RealArray& pos = initial_positions[pi];

                // convert position to index
                for (int dim=0; dim<AMREX_SPACEDIM; ++dim) {
                    pos_idx[dim] = std::floor((pos[dim] - plo[dim])*dxi[dim]);
                }

                if (box.contains(pos_idx)) {
                    Array<Real,3>& vel = initial_velocities[pi];
                    CParticleType p;
                    p.id()   = CParticleType::NextID();
                    p.cpu()  = ParallelDescriptor::MyProc();
                    AMREX_D_TERM(
                                p.pos(0) = pos[0];,
                            p.pos(1) = pos[1];,
                    p.pos(2) = pos[2];
                    )

                    p.rdata(+ParticleIdxR::VX) = vel[0];
                    p.rdata(+ParticleIdxR::VY) = vel[1];
                    p.rdata(+ParticleIdxR::VZ) = vel[2];

                    p.rdata(+ParticleIdxR::Charge) = charge;
                    p.rdata(+ParticleIdxR::Mass) = mass;

                    pc.push_back(p);

                    done[pi] = 1;
                }
            }
        }

        particles->Redistribute();
    }
}

void ChargedParticle::checkpoint(const std::string& dir)
{
    particles->Checkpoint(dir, "Particles_"+name, true, particle_real_names, particle_int_names);
}

void ChargedParticle::restart(const std::string& dir)
{
    particles->Restart(dir, "Particles_"+name);
}

void ChargedParticle::redistribute(int level, int finest_level, int ngrow)
{
    particles->Redistribute(level, finest_level, ngrow);
}

void ChargedParticle::clear()
{
    particles.reset();
}

void ChargedParticle::push_particles(MFIter& mfi,
                                     const FArrayBox& prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                     const int E_idx,
                                     const int B_idx,
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

    CParTileType& ptile = particles->ParticlesAt(level, mfi);

//    const Real          strttime = amrex::second();
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();


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

    Real scale_factor = 0.5*dt/MFP::Larmor;

    amrex::ParallelFor(n,
                       [=] AMREX_GPU_DEVICE (int i)
    {
        CParticleType& p  = p_pbox[i];
        if (p.id() <= 0) return;
        Real E[3] = {0.0,0.0,0.0};
        Real B[3] = {0.0,0.0,0.0};

        // implement particle boundary conditions here??

//        int valid = 0;

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

        // get the E and B fields at the particle according to the local slopes
        // obtained via the reconstructed cell face values

        for (int d1=0; d1<3; ++d1) {
            E[d1] = p4(iloc[0], iloc[1], iloc[2], E_idx + d1); // cell centre value
            B[d1] = p4(iloc[0], iloc[1], iloc[2], B_idx + d1); // cell centre value
            for (int d2=0; d2<AMREX_SPACEDIM; ++d2) {
                // adjustment due to slopes
                const Real E_lo_val = rlo4[d2](iloc[0], iloc[1], iloc[2],E_idx + d1);
                const Real E_hi_val = rhi4[d2](iloc[0], iloc[1], iloc[2],E_idx + d1);
                E[d1] += loc[d2]*(E_hi_val - E_lo_val);

                const Real B_lo_val = rlo4[d2](iloc[0], iloc[1], iloc[2],B_idx + d1);
                const Real B_hi_val = rhi4[d2](iloc[0], iloc[1], iloc[2],B_idx + d1);
                B[d1] += loc[d2]*(B_hi_val - B_lo_val);
            }
        }

        // apply updates to particle position and velocity using the Boris pusher

        Real q = p.rdata(+ParticleIdxR::Charge);
        Real m = p.rdata(+ParticleIdxR::Mass);

        Real E0[3], B0[3], Vn[3], V0[3], V1_[3], V1[3];
        for (int i=0; i<3; ++i) {
            E0[i] = scale_factor*q*E[i]/m;
            B0[i] = scale_factor*q*B[i]/m;
            Vn[i] = p.rdata(+ParticleIdxR::VX+i);
            V0[i] = Vn[i] + E0[i];
        }

        Real B02 = 1 + B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2];

        V1_[0] = 2*(V0[0] + V0[1]*B0[2] - V0[2]*B0[1])/B02;
        V1_[1] = 2*(V0[1] + V0[2]*B0[0] - V0[0]*B0[2])/B02;
        V1_[2] = 2*(V0[2] + V0[0]*B0[1] - V0[1]*B0[0])/B02;

        V1[0] = V0[0] + V1_[1]*B0[2] - V1_[2]*B0[1];
        V1[1] = V0[1] + V1_[2]*B0[0] - V1_[0]*B0[2];
        V1[2] = V0[2] + V1_[0]*B0[1] - V1_[1]*B0[0];

        for (int i=0; i<3; ++i) {
            p.rdata(+ParticleIdxR::VX+i) = V1[i] + E0[i];
        }

        for (int dim=0; dim<AMREX_SPACEDIM; ++dim) {
            p.pos(dim) += dt*p.rdata(+ParticleIdxR::VX+dim);
        }

    });

    return;
}

void ChargedParticle::calculate_source(MFIter& mfi, FArrayBox& S, Geometry& geom, int level) const
{

    const Real* dom_lo = geom.ProbLo();
    const Real* dxi = geom.InvCellSize();

    // zero out first
    S.setVal(0.0);

    const Real factor = MFP::Larmor/(MFP::lightspeed*MFP::Debye*MFP::Debye);

    Array4<Real> const& S4 = S.array();

    CParTileType& pc = particles->DefineAndReturnParticleTile(level, mfi.index(), mfi.LocalTileIndex());
    CParticleType *  AMREX_RESTRICT particle = &(pc.GetArrayOfStructs()[0]);
    const int np = pc.numParticles();

    for (int i=0; i<np; ++i) {

        const Real q = particle[i].rdata(+ParticleIdxR::Charge);

        Array<Real,3> vel;
        vel[0] = particle[i].rdata(+ParticleIdxR::VX);
        vel[1] = particle[i].rdata(+ParticleIdxR::VY);
        vel[2] = particle[i].rdata(+ParticleIdxR::VZ);

        // get the local coordinates of the particle
        Array<int,3> coord = {0,0,0};
        AMREX_D_TERM(
                    coord[0]=floor((particle[i].pos(0) - dom_lo[0])*dxi[0]);,
        coord[1]=floor((particle[i].pos(1) - dom_lo[1])*dxi[1]);,
        coord[2]=floor((particle[i].pos(2) - dom_lo[2])*dxi[2]);
        )

        // calculate the current
        for (int vi=0; vi<3; ++vi) {
            S4(coord[0], coord[1], coord[2], vi) -= factor*q*vel[vi];
        }
    }

    return;
}

#ifdef AMREX_USE_EB
void ChargedParticle::set_eb_bc(const sol::table &bc_def)
{

    std::string bc_type = bc_def.get<std::string>("type");


}
#endif

void ChargedParticle::write_info(nlohmann::json& js) const
{
    LagrangianState::write_info(js);

    js["state_idx"] = state_idx;
    js["type"] = tag;

}

#endif
