#ifndef MFP_HYDRO_BC_H
#define MFP_HYDRO_BC_H
#ifdef AMREX_USE_EB

#include "sol.hpp"
#include "map"
#include "MFP_hydro_viscous.H"

class HydroBoundaryEB
{
public:
    HydroBoundaryEB();
    ~HydroBoundaryEB();

    enum EBType {
        Dirichlet = 0,
        Slip,
        NoSlip
    };

    virtual const EBType get_type() const = 0;

    HydroRiemannSolver* flux_solver = nullptr;

    EBType type;
};

// custom boundary state where the user is able to optionally define any component
// of the primitive state vector
// Note that x points into the wall

class DirichletWall : public HydroBoundaryEB
{
public:
    DirichletWall();
    DirichletWall(HydroRiemannSolver* flux,
                  const Vector<std::string> &names,
                  const Vector<int> &vec_idx,
                  const sol::table &bc_def);
    ~DirichletWall();

    void solve(Array<Array<Real,3>,3> &wall_coord,
                      Array<Real,+HydroDef::FluxIdx::NUM> &state,
                      Array<Array<Real,+HydroDef::ConsIdx::NUM>, AMREX_SPACEDIM> &F,
                      const Real *dx);


    virtual const EBType get_type() const override {return EBType::Dirichlet;}

    Vector<std::pair<int, Real>> wall_value;

    static std::string tag;
};

/*
//-----------------------------------------------------------------------------

// standard slip wall

class HydroSlipWall : public HydroBoundaryEB
{
public:
    HydroSlipWall();
    HydroSlipWall(HydroRiemannSolver* flux);
    ~HydroSlipWall();

    void solve(Array<Array<Real,3>,3> &wall_coord,
                      Array<Real,+HydroDef::ConsIdx::NUM> &state,
                      Array<Array<Real,+HydroDef::ConsIdx::NUM>, AMREX_SPACEDIM> &F,
                      const Real *dx);

    virtual const EBType get_type() const override {return type;}

    EBType type = EBType::Slip;
    static std::string tag;

};

//-----------------------------------------------------------------------------

// no slip wall using the diffusive fluxes defined by the diffusion model used
// by the hydro state that it is interacting with

class HydroNoSlipWall : public HydroBoundaryEB
{
public:
    HydroNoSlipWall();
    HydroNoSlipWall(HydroRiemannSolver* flux, HydroViscous* visc, const sol::table &bc_def);
    ~HydroNoSlipWall();

    void solve(Array<Array<Real,3>,3> &wall_coord,
                      Array<Real, +HydroDef::ConsIdx::NUM> &state,
                      Array<Real,3> &normal_slope,
                      Array<Vector<Real>, AMREX_SPACEDIM> &F,
                      const Real* dx);

    virtual const EBType get_type() const override {return EBType::NoSlip;}

    HydroViscous* viscous = nullptr;
    Array<int,4> slopes_idx;
    Vector<Real> wall_velocity;
    Real wall_temp;

    static std::string tag;
};
*/
#endif
#endif // MFP_HYDRO_BC_H
