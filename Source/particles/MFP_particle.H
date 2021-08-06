#ifndef PARTICLE_H
#define PARTICLE_H
#ifdef AMREX_PARTICLES

#include <AMReX_AmrParticles.H>

#include "json.hpp"
#include "MFP_factory.H"
#include "MFP_eb_bc.H"

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
    virtual int num_source() {return 0;}
    virtual void calculate_source(MFIter& mfi, FArrayBox& S, Geometry &geom, int level) const {return;}
    virtual void clear() = 0;

    std::string name;
    int global_idx;
    Vector<std::pair<int,int>> associated_sources; // (system idx, src idx)

#ifdef AMREX_USE_EB
        bool eb_all_regular;
        int eb_index;
        const EB2::IndexSpace* eb2_index;
        Vector<std::unique_ptr<BoundaryEB>> eb_bcs;

        virtual void set_eb_bc(const sol::table &bc_def);
#endif
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
