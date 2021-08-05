#ifndef CHARGEDPARTICLE_H
#define CHARGEDPARTICLE_H
#ifdef AMREX_PARTICLES

#include "MFP_particle.H"
#include "MFP_utility.H"
#include "MFP_defs.H"

using namespace amrex;

class ChargedParticle : public ParticleMFP
{
public:
    ChargedParticle();
    ChargedParticle(const sol::table& def);
    virtual void init(AmrCore *amr_core, bool make_particles=true) override;

    virtual void checkpoint(const std::string& dir) override;
    virtual void restart(const std::string& dir) override;
    virtual void redistribute(int level, int finest_level, int ngrow) override;
    virtual void clear() override;

    void push_particles(MFIter& mfi,
                        const FArrayBox& prim,
                        Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                        Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                        const int E_idx,
                        const int B_idx,
                        const Real dt,
                        const Geometry geom,
                        const int level
                        EB_OPTIONAL(,const EBCellFlagFab& flag));

    virtual ParticleType get_type() const override {return ParticleMFP::ParticleType::Charged;}
    virtual const std::string& get_tag() const override {return tag;}
    virtual void write_info(nlohmann::json& js) const override;

    enum class ParticleIdxR : int {
        Mass=0,
        Charge,
        VX,
        VY,
        VZ,
        NUM
    };

    enum class ParticleIdxI : int {
        NUM=0
    };

    static std::string tag;
    static bool registered;

    int state_idx;
    int verbosity;
    Real charge, mass;
    Vector<RealArray> initial_positions;
    Vector<Array<Real,3>> initial_velocities;

    static Vector<std::string> particle_real_names;
    static Vector<std::string> particle_int_names;

    std::unique_ptr<AmrParticleContainer<+ParticleIdxR::NUM, +ParticleIdxI::NUM>> particles;
};

typedef AmrParticleContainer<+ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM> AmrCParContType;
typedef ParIter< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM, 0, 0 > CParIterType;
typedef ParticleContainer< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM, 0, 0 > CParContainerType;
typedef Particle< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM > CParticleType;
typedef ArrayOfStructs< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM > CParAoSType;
typedef ParticleTile< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM, 0, 0 > CParTileType;

#endif
#endif // CHARGEDPARTICLE_H
