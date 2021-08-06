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

#ifdef AMREX_USE_EB

void ParticleMFP::set_eb_bc(const sol::table &bc_def)
{
    return;
}
#endif

void ParticleMFP::write_info(nlohmann::json& js) const
{
    js["name"] = name;
    js["global_idx"] = global_idx;
}
#endif
