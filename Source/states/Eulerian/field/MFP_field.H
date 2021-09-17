#ifndef FIELDSTATE_H
#define FIELDSTATE_H

#include "MFP_eulerian.H"
#include "MFP_field_defs.H"
#include "MFP_state.H"
#include "MFP.H"
#include "MFP_bc.H"
#include "MFP_field_bc.H"
#include "MFP_field_riemann.H"

class FieldState : public EulerianState
{
public:
    FieldState();
    FieldState(const sol::table &def);
    ~FieldState();

    static bool registered;
    static std::string tag;


    static std::map<std::string, int> bc_names;

    static Vector<set_bc> bc_set;

    static Vector<std::string> cons_names;
    static Array<int,+FieldDef::VectorIdx::Cons> vector_idx;

    Real div_speed;
    Real relative_div_speed;
    Real fastest_speed;
    bool is_static;
    int project_divergence;

    std::map<std::string,Optional3D1VFunction> functions;

    std::unique_ptr<FieldRiemannSolver> flux_solver;

    virtual StateType get_type() const override {return StateType::Field;}
    virtual int get_num_cons() const override {return +FieldDef::ConsIdx::NUM;}
    virtual const Vector<std::string>& get_cons_names() const override {return cons_names;}

    virtual Vector<std::string> get_plot_output_names() const override {return cons_names;}
    virtual void get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                             #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                             #endif
                                 ) const override;

    virtual void variable_setup() override;
    virtual void init_data(MFP *mfp) override;

    void set_udf();

    virtual void init_from_lua() override;

#ifdef AMREX_USE_EB
    virtual void set_eb_bc(const sol::table &bc_def) override;
    Vector<std::unique_ptr<FieldBoundaryEB>> eb_bcs;
#endif

    virtual Real get_allowed_time_step(MFP* mfp) const override;

    static FieldState& get_state_global(const int idx) {return static_cast<FieldState&>(MFP::get_state(idx));}

};

#endif // FIELDSTATE_H
