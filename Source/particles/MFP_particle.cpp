#ifdef AMREX_PARTICLES
#include "MFP_particle.H"

ParticleMFP::ParticleMFP()
{
}

ParticleMFP::~ParticleMFP()
{
}

ClassFactory<ParticleMFP>& GetParticleFactory()
{
    static ClassFactory<ParticleMFP> F;
    return F;
}

void ParticleMFP::write_info(nlohmann::json& js) const
{
    js["name"] = name;
    js["global_idx"] = global_idx;
}
#endif
