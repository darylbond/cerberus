#ifndef TRACER_H
#define TRACER_H
#ifdef AMREX_PARTICLES

#include "MFP_particle.H"
#include "MFP_global.H"

using namespace amrex;

class TracerParticle : public ParticleMFP
{
public:
    TracerParticle();
    TracerParticle(const sol::table& def);
    virtual void init(AmrCore *amr_core, bool make_particles=true) override;

    virtual void checkpoint(const std::string& dir) override;
    virtual void restart(const std::string& dir) override;
    virtual void redistribute(int level, int finest_level, int ngrow) override;
    virtual void clear() override;

    void push_particles(MFIter& mfi,
                         const FArrayBox& prim,
                         Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                         Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                         const int vel_idx,
                         const Real dt,
                         const Geometry geom,
                        const int level
                         EB_OPTIONAL(,const EBCellFlagFab& flag));

    virtual ParticleType get_type() const override {return ParticleMFP::ParticleType::Tracer;}
    virtual const std::string& get_tag() const override {return tag;}
    virtual void write_info(nlohmann::json& js) const override;

    enum class ParticleIdxR : int {
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

    static Vector<std::string> particle_real_names;
    static Vector<std::string> particle_int_names;

    Vector<RealArray> initial_positions;
    int verbosity;



    std::unique_ptr<AmrParticleContainer<+ParticleIdxR::NUM, +ParticleIdxI::NUM>> particles;
};

typedef AmrParticleContainer<+TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM> AmrTParContType;
typedef ParIter< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM, 0, 0 > TParIterType;
typedef ParticleContainer< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM, 0, 0 > TParContainerType;
typedef Particle< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM > TParticleType;
typedef ArrayOfStructs< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM > TParAoSType;
typedef ParticleTile< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM, 0, 0 > TParTileType;

#endif
#endif // TRACER_H
