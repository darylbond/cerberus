#include "MFP_plasma5.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"
#include "Eigen"
#include "Dense"

std::string Plasma5::tag = "plasma5";
bool Plasma5::registered = GetSourceFactory().Register(Plasma5::tag, SourceBuilder<Plasma5>);
Plasma5::TimeIntegrator Plasma5::solver;

Plasma5::Plasma5(){}
Plasma5::~Plasma5(){}

Plasma5::Plasma5(const int idx, const sol::table &def)
{
    src_idx = idx;
    name = def["name"];

    std::string solver_name = def["solver"];

    if (solver_name == "explicit") {
        solver = TimeIntegrator::ForwardsEuler;
    } else if (solver_name == "implicit") {
        solver = TimeIntegrator::BackwardsEuler;
    } else {
        Abort("Plasma5 source needs to know what type of solver to use, options are ['explicit','implicit']");
    }

    const sol::table states = def["states"];

    for (const auto& key_value_pair : states) {
        std::string state_name = key_value_pair.second.as<std::string>();
        State& istate = MFP::get_state(state_name);

        switch (istate.get_type()) {
        case State::StateType::Field:
            if (field != nullptr) Abort("Only one field state can be set for the Plasma5 source "+name);
            field = static_cast<FieldState*>(&istate);
            break;
        case State::StateType::Hydro:
            species.push_back(static_cast<HydroState*>(&istate));
            break;
        default:
            Abort("An invalid state has been defined for the Plasma5 source "+name);
        }
    }
    return;
}

void Plasma5::solve(MFP* mfp, const Real dt) const
{
    BL_PROFILE("Plasma5::solve");

    switch(solver) {
    case TimeIntegrator::ForwardsEuler :
        explicit_solve(mfp, dt);
        break;
    case TimeIntegrator::BackwardsEuler :
        implicit_solve(mfp,dt);
        break;
    default:
        Abort("How did we get here?");
    }


}


void Plasma5::explicit_solve(MFP* mfp, const Real dt) const
{
    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    MultiFab* field_data = &(mfp->get_new_data(field->data_idx));


    size_t n_species = species.size();

    Vector<MultiFab*> species_data;
    for (const HydroState* hstate : species) {
        species_data.push_back(&(mfp->get_new_data(hstate->data_idx)));
    }

    Vector<Array4<Real>> species4(n_species);

    // define some 'registers'

    Real pD, pB;
    Real Bx, By, Bz;
    Real Ex, Ey, Ez;
    Real ep;

    Real q, m, r;
    Real rho, mx, my, mz, nrg, alpha;
    Real u, v, w;

    const Real Larmor = MFP::Larmor;
    const Real Debye = MFP::Debye;
    const Real lightspeed = MFP::lightspeed;

    const Real c_h = field->div_speed;
    const Real div_damp = field->div_damping;
    const Real c_d = div_damp != 0.0 ? c_h/div_damp : 0.0;

    const Real f1 = dt*Larmor/(Debye*Debye*lightspeed);
    const Real f2 = dt*c_h*c_h*f1/lightspeed;
    const Real f3 = c_d != 0.0 ? dt*c_h*c_h/(c_d*c_d) : 0.0;

    // get charge and current density
    Real charge_density, current_x, current_y, current_z;



    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(field->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        //        bool skip = false;
        //        for (const HydroState* hstate : species) {
        //            EBData& eb = mfp->get_eb_data(hstate->global_idx);
        //            const FArrayBox& vfrac_hydro = eb.volfrac[mfi];
        //            if (vfrac_hydro.getType() == FabType::covered) {
        //                skip = true;
        //                break;
        //            }
        //            //            if (???) {
        //            //                Abort("All EB data is not the same for source "+name);
        //            //            }
        //        }
        //        if (skip) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        Array4<Real> const& field4 = field_data->array(mfi);

        for (int n=0; n<n_species; ++n) {
            species4[n] = species_data[n]->array(mfi);
        }


        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif


                    ep = field4(i,j,k,+FieldDef::ConsIdx::ep);

                    // magnetic field
                    Bx = field4(i,j,k,+FieldDef::ConsIdx::Bx);
                    By = field4(i,j,k,+FieldDef::ConsIdx::By);
                    Bz = field4(i,j,k,+FieldDef::ConsIdx::Bz);
                    pB = field4(i,j,k,+FieldDef::ConsIdx::psi);

                    // electric field
                    Ex = field4(i,j,k,+FieldDef::ConsIdx::Dx)/ep;
                    Ey = field4(i,j,k,+FieldDef::ConsIdx::Dy)/ep;
                    Ez = field4(i,j,k,+FieldDef::ConsIdx::Dz)/ep;
                    pD = field4(i,j,k,+FieldDef::ConsIdx::phi);

                    // get charge and current density
                    charge_density = 0.0;
                    current_x = 0.0;
                    current_y = 0.0;
                    current_z = 0.0;

                    for (size_t n = 0; n < n_species; ++n) {

                        rho =   species4[n](i,j,k,+HydroDef::ConsIdx::Density);
                        mx =    species4[n](i,j,k,+HydroDef::ConsIdx::Xmom);
                        my =    species4[n](i,j,k,+HydroDef::ConsIdx::Ymom);
                        mz =    species4[n](i,j,k,+HydroDef::ConsIdx::Zmom);
                        nrg =   species4[n](i,j,k,+HydroDef::ConsIdx::Eden);
                        alpha = species4[n](i,j,k,+HydroDef::ConsIdx::Tracer)/rho;

                        m = species[n]->get_mass(alpha);
                        q = species[n]->get_charge(alpha);

                        r = q/m;

                        charge_density += rho*r;
                        current_x += r*mx;
                        current_y += r*my;
                        current_z += r*mz;

                        u = mx/rho;
                        v = my/rho;
                        w = mz/rho;

                        species4[n](i,j,k,+HydroDef::ConsIdx::Xmom) += dt*(rho*r/Larmor)*(lightspeed*Ex + v*Bz - w*By);
                        species4[n](i,j,k,+HydroDef::ConsIdx::Ymom) += dt*(rho*r/Larmor)*(lightspeed*Ey + w*Bx - u*Bz);
                        species4[n](i,j,k,+HydroDef::ConsIdx::Zmom) += dt*(rho*r/Larmor)*(lightspeed*Ez + u*By - v*Bx);
                        species4[n](i,j,k,+HydroDef::ConsIdx::Eden) += dt*(rho*r*lightspeed)/Larmor*(u*Ex + v*Ey + w*Ez);
                    }

                    // electric field and divergence constraint sources



                    field4(i,j,k,+FieldDef::ConsIdx::Dx)  -= f1*current_x;
                    field4(i,j,k,+FieldDef::ConsIdx::Dy)  -= f1*current_y;
                    field4(i,j,k,+FieldDef::ConsIdx::Dz)  -= f1*current_z;

                    // source for D divergence correction
                    field4(i,j,k,+FieldDef::ConsIdx::phi) *= (1-f3);
                    field4(i,j,k,+FieldDef::ConsIdx::phi) +=  f2*charge_density;

                    // source for B divergence correction
                    field4(i,j,k,+FieldDef::ConsIdx::psi) *= (1 - f3);
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}

void Plasma5::implicit_solve(MFP* mfp, const Real dt) const
{

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    MultiFab* field_data = &(mfp->get_new_data(field->data_idx));


    size_t n_species = species.size();

    Vector<MultiFab*> species_data;
    for (const HydroState* hstate : species) {
        species_data.push_back(&(mfp->get_new_data(hstate->data_idx)));
    }

    Vector<Array4<Real>> species4(n_species);


    // define the linear system matrix/vector and the associated solver
    int N = 3*(n_species+1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(N,N);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd c(N);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> cphqr;

    Real trace_B, trace_C, mc;

    Real pD, pB;
    Real Bx, By, Bz;
    Real Ex, Ey, Ez;
    Real ep;

    Real charge_density;
    Real q, m, r, g;
    Real rho, mx, my, mz, alpha;
    Vector<Real> R(n_species);

    const Real Larmor = MFP::Larmor;
    const Real Debye = MFP::Debye;
    const Real lightspeed = MFP::lightspeed;
    const Real c_h = field->div_speed;
    const Real div_damp = field->div_damping;
    const Real c_d = div_damp != 0.0 ? c_h/div_damp : 0.0;

    const Real f1 = dt*Larmor/(Debye*Debye*lightspeed);
    const Real f2 = dt*c_h*c_h*f1/lightspeed;
    const Real f3 = c_d != 0.0 ? dt*c_h*c_h/(c_d*c_d) : 0.0;


    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(field->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        Array4<Real> const& field4 = field_data->array(mfi);

        for (int n=0; n<n_species; ++n) {
            species4[n] = species_data[n]->array(mfi);
        }


        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif


                    ep = field4(i,j,k,+FieldDef::ConsIdx::ep);

                    // magnetic field
                    Bx = field4(i,j,k,+FieldDef::ConsIdx::Bx);
                    By = field4(i,j,k,+FieldDef::ConsIdx::By);
                    Bz = field4(i,j,k,+FieldDef::ConsIdx::Bz);
                    pB = field4(i,j,k,+FieldDef::ConsIdx::psi);

                    // electric field
                    Ex = field4(i,j,k,+FieldDef::ConsIdx::Dx)/ep;
                    Ey = field4(i,j,k,+FieldDef::ConsIdx::Dy)/ep;
                    Ez = field4(i,j,k,+FieldDef::ConsIdx::Dz)/ep;
                    pD = field4(i,j,k,+FieldDef::ConsIdx::phi);


                    // load in electric field
                    b(n_species*3 + 0) = Ex;
                    b(n_species*3 + 1) = Ey;
                    b(n_species*3 + 2) = Ez;

                    charge_density = 0.0;

                    for (size_t n = 0; n < n_species; ++n) {

                        rho =   species4[n](i,j,k,+HydroDef::ConsIdx::Density);
                        mx =    species4[n](i,j,k,+HydroDef::ConsIdx::Xmom);
                        my =    species4[n](i,j,k,+HydroDef::ConsIdx::Ymom);
                        mz =    species4[n](i,j,k,+HydroDef::ConsIdx::Zmom);
                        alpha = species4[n](i,j,k,+HydroDef::ConsIdx::Tracer)/rho;

                        m = species[n]->get_mass(alpha);
                        q = species[n]->get_charge(alpha);

                        r = q/m;
                        R[n] = r;

                        charge_density += rho*r;

                        charge_density += rho*r;
                        trace_B = -dt*r*rho*lightspeed/Larmor;
                        trace_C = dt*r*Larmor/(lightspeed*Debye*Debye);
                        mc = dt*r/Larmor;

                        // diagonal elements
                        A(3*n + 0, 3*n + 1) = -Bz*mc;
                        A(3*n + 0, 3*n + 2) =  By*mc;

                        A(3*n + 1, 3*n + 0) =  Bz*mc;
                        A(3*n + 1, 3*n + 2) = -Bx*mc;

                        A(3*n + 2, 3*n + 0) = -By*mc;
                        A(3*n + 2, 3*n + 1) =  Bx*mc;

                        // end column elements
                        A(3*n + 0, 3*n_species + 0) = trace_B;
                        A(3*n + 1, 3*n_species + 1) = trace_B;
                        A(3*n + 2, 3*n_species + 2) = trace_B;

                        // end row elements
                        A(3*n_species + 0, 3*n + 0) = trace_C;
                        A(3*n_species + 1, 3*n + 1) = trace_C;
                        A(3*n_species + 2, 3*n + 2) = trace_C;

                        // b vector (momentum components)
                        b(n*3 + 0) = mx;
                        b(n*3 + 1) = my;
                        b(n*3 + 2) = mz;

                    }

                    // solve linear system for updated momentum and electric field
                    cphqr.compute(A);
                    c = cphqr.solve(b);


                    // the updated D field
                    field4(i,j,k,+FieldDef::ConsIdx::Dx) = c(n_species*3 + 0)*ep;
                    field4(i,j,k,+FieldDef::ConsIdx::Dy) = c(n_species*3 + 1)*ep;
                    field4(i,j,k,+FieldDef::ConsIdx::Dz) = c(n_species*3 + 2)*ep;

                    // new electric field
                    Ex = c(n_species*3 + 0);
                    Ey = c(n_species*3 + 1);
                    Ez = c(n_species*3 + 2);

                    // get the updates for hydro species
                    for (int n=0; n<n_species; ++n) {
                        species4[n](i,j,k,+HydroDef::ConsIdx::Xmom) = mx = c(n*3 + 0);
                        species4[n](i,j,k,+HydroDef::ConsIdx::Ymom) = my = c(n*3 + 1);
                        species4[n](i,j,k,+HydroDef::ConsIdx::Zmom) = mz = c(n*3 + 2);

                        species4[n](i,j,k,+HydroDef::ConsIdx::Eden) += dt*(R[n]*lightspeed/Larmor)*(Ex*mx + Ey*my + Ez*mz);

                    }


                    // source for E divergence correction
                    field4(i,j,k,+FieldDef::ConsIdx::phi) *= (1-f3);
                    field4(i,j,k,+FieldDef::ConsIdx::phi) += f2*charge_density;

                    // source for B divergence correction
                    field4(i,j,k,+FieldDef::ConsIdx::psi) *= (1 - f3);

                }
            }
        }
        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
