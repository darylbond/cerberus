#ifndef MFP_PARTICLE_CURRENT_H
#define MFP_PARTICLE_CURRENT_H

#include "MFP_source.H"

#include "sol.hpp"

class ParticleCurrentSource : public SourceTerm
{
    public:
        ParticleCurrentSource();
        ParticleCurrentSource(const sol::table& def);
        ~ParticleCurrentSource();

        virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;

        static bool valid_state(const int global_idx);
        static bool valid_solver(const int solve_idx);


        virtual SourceType get_type() const override {return SourceType::ParticleCurrent;}
        virtual const std::string& get_tag() const override {return tag;}

        static std::string tag;
        static bool registered;
};

#endif // MFP_PARTICLE_CURRENT_H
