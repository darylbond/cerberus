#ifndef MFP_STATE_H
#define MFP_STATE_H

#include "forward.hpp"
#include "json.hpp"

#include "MFP_optional_func.H"
#include "MFP_utility.H"
#include "MFP_factory.H"
#include "MFP_refine.H"

class State;
class MFP;

#include "MFP_eb_structs.H"
#include "MFP_fillbc.H"


using namespace amrex;

class State
{
public:
    State();
    virtual ~State();

    std::string name;
    int global_idx;
    int data_idx;

    enum class StateClassification {
        Null = -1,
        Eulerian,
        Lagrangian,
        NUM
    };

    enum class StateType {
        Null = -1,
        Hydro,
        Field,
        TracerParticle,
        ChargedParticle,
        NUM
    };

    Vector<int> associated_sources;

#ifdef AMREX_USE_EB
    bool eb_all_regular;
    int eb_index;
    const EB2::IndexSpace* eb2_index;
#endif

    BoundaryState boundary_conditions;


    std::unique_ptr<Refinement> refine;
    void get_refinement_tags(MFP* mfp, TagBoxArray& tags);

    virtual StateClassification get_classification() const {return StateClassification::Null;}
    virtual StateType get_type() const {return StateType::Null;}

    virtual Vector<std::string> get_plot_output_names() const {return {};}
    virtual void get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                             #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                             #endif
                                 ) const {return;}

    virtual void init_from_lua() {}
    virtual void post_init_from_lua() {}
    virtual void variable_setup(Vector<int> periodic) {}
    virtual void init_data(MFP* mfp) {}

    virtual Real get_allowed_time_step(MFP* mfp) const {return 0.0;}


    Vector<BoundaryInfo> get_bc_limits(const Box& box,
                                       const Geometry &geom) const;

#ifdef AMREX_USE_EB

    virtual void set_eb_bc(const sol::table &bc_def) {}

    void update_eb_vfrac(const Geometry &geom,
                         FArrayBox& vfrac) const;

    void update_eb_flags(const Geometry &geom,
                         EBCellFlagFab& flags) const;

    // keep track of the embedded boundaries that are being applied by this state
    // first index is the state local definition of the eb bc
    // second index is for grabbing any data to go with the bc
    Vector<std::pair<int,int>> eb_bc_index;

#endif

    virtual void write_info(nlohmann::json& js) const;

    virtual std::string str() const;

};

template <typename D>
std::unique_ptr<State> StateBuilder(const sol::table& def)
{
    std::string type = def["type"];
    if (type == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

ClassFactory<State>& GetStateFactory();

#endif // MFP_GENERIC_H
