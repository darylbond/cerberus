#ifndef HYDROTRACER_H
#define HYDROTRACER_H

#include "MFP_source.H"
#include "MFP_hydro.H"
#include "MFP_tracer.H"

class HydroTracer : public Source
{
public:
    HydroTracer();
    HydroTracer(const int idx, const sol::table &def);
    ~HydroTracer();

    void push_particles(TParTileType& ptile,
                        const FArrayBox& prim,
                        const Geometry geom,
                        const Real dt
                        #ifdef AMREX_USE_EB
                        ,const EBCellFlagFab& flag
                        #endif
                        ) const;

    virtual void solve(MFP* mfp, const Real dt) const override;


    virtual std::string get_tag() const override {return tag;}

    HydroState* hydro_state;
    TracerParticle* tracer_state;

    static bool registered;
    static std::string tag;
};

#endif // TRACER_H
