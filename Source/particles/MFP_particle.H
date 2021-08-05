#ifndef PARTICLE_H
#define PARTICLE_H
#ifdef AMREX_PARTICLES

#include <AMReX_AmrParticles.H>

#include "json.hpp"
#include "MFP_factory.H"

using namespace amrex;


class ParticleMFP
{
public:
    ParticleMFP();
    virtual ~ParticleMFP();

    enum class ParticleType : int {
        Tracer,
        Charged
    };

    virtual const std::string& get_tag() const = 0;
    virtual ParticleType get_type() const = 0;
    virtual void write_info(nlohmann::json& js) const;
    virtual void init(AmrCore *amr_core, bool make_particles=true) = 0;
    virtual void checkpoint(const std::string& dir) = 0;
    virtual void restart(const std::string& dir) = 0;
    virtual void redistribute(int level, int finest_level, int ngrow) = 0;
    virtual void clear() = 0;

    std::string name;
    int global_idx;
};

template <typename D>
std::unique_ptr<ParticleMFP> ParticleBuilder(const sol::table& def)
{
    if (def["type"] == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

ClassFactory<ParticleMFP>& GetParticleFactory();

#endif
#endif // PARTICLE_H
