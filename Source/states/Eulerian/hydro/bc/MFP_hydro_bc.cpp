#ifdef AMREX_USE_EB

#include "MFP_hydro_bc.H"
#include "MFP_hydro.H"
#include "MFP_transforms.H"


//-----------------------------------------------------------------------------

std::string DirichletWall::tag = "dirichlet";

DirichletWall::DirichletWall(){}
DirichletWall::~DirichletWall(){}

DirichletWall::DirichletWall(HydroRiemannSolver *flux,
                             const Vector<std::string> &names,
                             const Vector<int> &vec_idx,
                             const sol::table &bc_def)
{
    BL_PROFILE("DirichletWall::DirichletWall");
    flux_solver = flux;

    // grab the wall state from the lua definition
    // only get stuff that is defined, anything else will
    // be filled in from the fluid side of the wall
    int i = 0;
    for (const auto &name : names) {
        sol::object val = bc_def[name];
        if (val.valid()) {
            wall_value.push_back({i, val.as<Real>()});
        }
        ++i;
    }
}

void DirichletWall::solve(Array<Array<Real,3>,3> &wall_coord,
                          Array<Real,+HydroDef::FluxIdx::NUM> &state,
                          Array<Array<Real,+HydroDef::ConsIdx::NUM>, AMREX_SPACEDIM> &F,
                          const Real *dx)
{
    BL_PROFILE("DirichletWall::solve");
    //
    // get the inviscid flux
    //

    transform_global2local(state, wall_coord, HydroState::flux_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    Array<Real, +HydroDef::FluxIdx::NUM> W = state;

    for (const auto& pair : wall_value) {
        W[pair.first] = pair.second;
    }

    Array<Real,+HydroDef::ConsIdx::NUM> normal_flux;
    flux_solver->solve(state, W, normal_flux, nullptr);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, HydroState::cons_vector_idx);

    // split it up into components
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for (int n=0; n<normal_flux.size(); ++n) {
            F[d][n] = wall_coord[0][d]*normal_flux[n];
        }
    }

}
/*
//-----------------------------------------------------------------------------

std::string HydroSlipWall::tag = "slip_wall";

HydroSlipWall::HydroSlipWall(){}
HydroSlipWall::~HydroSlipWall(){}

HydroSlipWall::HydroSlipWall(HydroRiemannSolver *flux)
{
    flux_solver = flux;
}

void HydroSlipWall::solve(Array<Array<Real,3>,3> &wall_coord,
                          Array<Real, +HydroDef::ConsIdx::NUM> &state,
                          Array<Array<Real, +HydroDef::ConsIdx::NUM>, AMREX_SPACEDIM> &F,
                          const Real* dx)
{
    BL_PROFILE("HydroSlipWall::solve");
    //
    // get the inviscid flux
    //

    transform_global2local(state, wall_coord, HydroState::cons_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    Vector<Real> W = state;

    W[+HydroDef::PrimIdx::Xvel] *= -1;

    Vector<Real> normal_flux(+HydroDef::ConsIdx::NUM);

    Real shk = 1.0; // consider there to be a shock present if we have a switching flux
    flux_solver->solve(state, W, normal_flux, &shk);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, HydroState::cons_vector_idx);

    // split it up into components
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
            F[d][n] = wall_coord[0][d]*normal_flux[n];
        }
    }

    return;
}

//-----------------------------------------------------------------------------

std::string HydroNoSlipWall::tag = "no_slip_wall";

HydroNoSlipWall::HydroNoSlipWall(){}
HydroNoSlipWall::~HydroNoSlipWall(){}

HydroNoSlipWall::HydroNoSlipWall(HydroRiemannSolver *flux,
                                 HydroViscous *visc,
                                 const sol::table &bc_def)
{
    BL_PROFILE("HydroNoSlipWall::HydroNoSlipWall");
    flux_solver = flux;
    viscous = visc;

    // grab the wall velocity from the lua definition
    wall_velocity.resize(3,0.0);
    wall_velocity[1] = bc_def["v1"].get_or(0.0);
    wall_velocity[2] = bc_def["v2"].get_or(0.0);

    // grab a wall temperature
    wall_temp = bc_def["T"].get_or(-1.0);

    slopes_idx = {+HydroDef::PrimIdx::Xvel, +HydroDef::PrimIdx::Yvel, +HydroDef::PrimIdx::Zvel, -1};

    if (wall_temp > 0.0) {
        slopes_idx[3] = +HydroDef::PrimIdx::Temp;
    }

}

void HydroNoSlipWall::solve(Array<Array<Real,3>,3> &wall_coord,
                            Array<Real,+HydroDef::ConsIdx::NUM> &state,
                            Array<Real,3> &normal_slope,
                            Array<Vector<Real>, AMREX_SPACEDIM> &F,
                            const Real *dx)
{
    BL_PROFILE("HydroNoSlipWall::solve");
    //
    // get the inviscid flux
    //

    transform_global2local(state, wall_coord,  HydroState::flux_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    Array<Real,+HydroDef::PrimIdx::NUM> W = state;

    W[+HydroDef::PrimIdx::Xvel] *= -1;

    Vector<Real> normal_flux(flux_solver->get_flux_size());

    Real shk = 1.0; // consider there to be a shock present if we have a switching flux
    flux_solver->solve(state, W, normal_flux, &shk);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord,  HydroState::cons_vector_idx);

    // split it up into components
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for (int n=0; n<normal_flux.size(); ++n) {
            F[d][n] = wall_coord[0][d]*normal_flux[n];
        }
    }

    //
    // get the viscous flux
    //

    if (viscous) {

        Array<Real, 3> &wall_normal = wall_coord[0];

        Real T, mu, kappa;
        viscous->get_neutral_coeffs(state, T, mu, kappa);

        Array<Real, 3> dxinv = {0,0,0};
        AMREX_D_TERM(dxinv[0] = 1/dx[0];,dxinv[1] = 1/dx[1];,dxinv[2] = 1/dx[2];)

        //
        // transform them to d/dx, d/dy and d/dz given transverse derivatives are zero
        Real dudx = normal_slope[0] * wall_normal[0]*dxinv[0];
        Real dudy = normal_slope[0] * wall_normal[1]*dxinv[1];
        Real dudz = normal_slope[0] * wall_normal[2]*dxinv[2];
        //
        Real dvdx = normal_slope[1] * wall_normal[0]*dxinv[0];
        Real dvdy = normal_slope[1] * wall_normal[1]*dxinv[1];
        Real dvdz = normal_slope[1] * wall_normal[2]*dxinv[2];
        //
        Real dwdx = normal_slope[2] * wall_normal[0]*dxinv[0];
        Real dwdy = normal_slope[2] * wall_normal[1]*dxinv[1];
        Real dwdz = normal_slope[2] * wall_normal[2]*dxinv[2];

        Real divu = dudx+dvdy+dwdz;

        Real tautmp = 2.0/3.0*mu*divu;

        Real tauxx = mu*2*dudx + tautmp;
        Real tauyy = mu*2*dvdy + tautmp;

        Real tauxy = mu*(dudy+dvdx);
        Real tauxz = mu*(dudz+dwdx);
        Real tauyz = mu*(dwdy+dvdz);

        Real xvel = wall_velocity[0]*wall_coord[0][0] + wall_velocity[1]*wall_coord[1][0] + wall_velocity[2]*wall_coord[2][0];
        Real yvel = wall_velocity[0]*wall_coord[0][1] + wall_velocity[1]*wall_coord[1][1] + wall_velocity[2]*wall_coord[2][1];
        Real zvel = wall_velocity[0]*wall_coord[0][2] + wall_velocity[1]*wall_coord[1][2] + wall_velocity[2]*wall_coord[2][2];

        F[0][+HydroDef::ConsIdx::Xmom] -= tauxx;
        F[0][+HydroDef::ConsIdx::Ymom] -= tauxy;
        F[0][+HydroDef::ConsIdx::Zmom] -= tauxz;
        F[0][+HydroDef::ConsIdx::Eden] -= xvel*tauxx + yvel*tauxy + zvel*tauxz;

        F[1][+HydroDef::ConsIdx::Xmom] -= tauxy;
        F[1][+HydroDef::ConsIdx::Ymom] -= tauyy;
        F[1][+HydroDef::ConsIdx::Zmom] -= tauyz;
        F[1][+HydroDef::ConsIdx::Eden] -= xvel*tauxy + yvel*tauyy + zvel*tauyz;

#if AMREX_SPACEDIM == 3
        Real tauzz = mu*2*dwdz + tautmp;
        F[2][+HydroDef::ConsIdx::Xmom] -= tauxz;
        F[2][+HydroDef::ConsIdx::Ymom] -= tauyz;
        F[2][+HydroDef::ConsIdx::Zmom] -= tauzz;
        F[2][+HydroDef::ConsIdx::Eden] -= xvel*tauxz + yvel*tauyz + zvel*tauzz;
#endif

        // if we haven't specified a wall temperature
        // assume adiabatic (dTdn = 0)
        if (wall_temp > 0.0) {
            // calculate the temperature gradient
            Real dTdx = normal_slope[3] * wall_normal[0]*dxinv[0];
            Real dTdy = normal_slope[3] * wall_normal[1]*dxinv[1];

            F[0][+HydroDef::ConsIdx::Eden] -= kappa*dTdx;
            F[1][+HydroDef::ConsIdx::Eden] -= kappa*dTdy;

#if AMREX_SPACEDIM == 3
            Real dTdz = normal_slope[3] * wall_normal[2]*dxinv[2];
            F[2][+HydroDef::ConsIdx::Eden] -= kappa*dTdz;
#endif
        }

    }

}
*/
#endif
