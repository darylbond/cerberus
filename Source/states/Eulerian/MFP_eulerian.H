#ifndef EULERIANSTATE_H
#define EULERIANSTATE_H

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_reconstruction.H"

class EulerianState : public State
{
public:
    EulerianState();
    ~EulerianState();

    int num_grow;
    bool reflux;

    std::map<std::string,Optional3D1VFunction> functions;
    std::map<int,Optional3D1VFunction*> dynamic_functions; // index is for primitive values

    virtual StateClassification get_classification() const override {return StateClassification::Eulerian;}
    virtual int get_num_cons() const {return 0;}
    virtual const Vector<std::string>& get_cons_names() const = 0;

    virtual void set_reconstruction();
    virtual void set_reflux();
    virtual void set_shock_detector();

#ifdef AMREX_USE_EB
    virtual void set_eb_divergence();
#endif

    virtual void init_from_lua() override;

    virtual bool is_transported() {return true;}

    virtual void calc_primitives(const Box& box,
                                 FArrayBox& cons,
                                 FArrayBox &prim,
                                 const Real* dx,
                                 const Real t,
                                 const Real* prob_lo
                             #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                             #endif
                                 ) const = 0;

    virtual void update_boundary_cells(const Box& box,
                                       const Geometry &geom,
                                       FArrayBox &prim,
                                   #ifdef AMREX_USE_EB
                                       const FArrayBox& vfrac,
                                   #endif
                                       const Real time) const = 0;

    static EulerianState& get_state(const int idx) {return static_cast<EulerianState&>(MFP::get_state(idx));}

};

#endif // EULERIANSTATE_H
