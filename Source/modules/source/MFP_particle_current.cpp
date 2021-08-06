#include "MFP_particle_current.H"

using GD = GlobalData;

//---------------------------------------------------------------------------------------------

std::string ParticleCurrentSource::tag = "particle_current";
bool ParticleCurrentSource::registered = GetSourceTermFactory().Register(ParticleCurrentSource::tag, SourceTermBuilder<ParticleCurrentSource>);

ParticleCurrentSource::ParticleCurrentSource(){}

ParticleCurrentSource::ParticleCurrentSource(const sol::table &def)
{

    name = def.get<std::string>("name");

    if (!ParticleCurrentSource::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &ParticleCurrentSource::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    return;
}

ParticleCurrentSource::~ParticleCurrentSource()
{
    // do nothing
}



int ParticleCurrentSource::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{
    BL_PROFILE("ParticleCurrentSource::fun_rhs");

    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        std::cout << "idx = [" << idx.global << ", " << idx.local << ", " << idx.solver << "]" << std::endl;


    }

    return 0;
}

bool ParticleCurrentSource::valid_state(const int global_idx)
{

    if (GD::idx_is_state(global_idx)) {
        State& istate = GD::get_state(global_idx);
        switch (istate.get_type()) {
            case +StateType::isField:
                return true;
            default:
                return false;
        }
    } else {
        ParticleMFP &ipart = GD::get_particles(global_idx);
        switch (ipart.get_type()) {
            case ParticleMFP::ParticleType::Charged:
                return true;
            default:
                return false;
        }
    }
}

bool ParticleCurrentSource::valid_solver(const int solve_idx)
{
    if (solve_idx != +SolveType::Explicit) {
        return false;
    }
    return true;
}
