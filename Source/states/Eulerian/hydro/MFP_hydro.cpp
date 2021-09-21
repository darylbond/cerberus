#include "MFP_hydro.H"
#include "MFP_lua.H"
#include "MFP.H"
#include "MFP_diagnostics.H"
#include "MFP_transforms.H"
#include "MFP_hydro_refine.H"

#include "Eigen"


Vector<std::string> HydroState::cons_names = {
    "rho",
    "x_mom",
    "y_mom",
    "z_mom",
    "nrg",
    "tracer",
    "shock",
};

Vector<std::string> HydroState::prim_names = {
    "rho",
    "x_vel",
    "y_vel",
    "z_vel",
    "p",
    "T",
    "alpha"
};

Array<int,1> HydroState::flux_vector_idx = {+HydroDef::FluxIdx::Xvel};
Array<int,1> HydroState::cons_vector_idx = {+HydroDef::ConsIdx::Xmom};
Array<int,1> HydroState::prim_vector_idx = {+HydroDef::PrimIdx::Xvel};

std::map<std::string, int> HydroState::bc_names = {{"interior",  PhysBCType::interior},
                                                   {"inflow",    PhysBCType::inflow},
                                                   {"outflow",   PhysBCType::outflow},
                                                   {"symmetry",  PhysBCType::symmetry},
                                                   {"slipwall",  PhysBCType::slipwall},
                                                   {"noslipwall",PhysBCType::noslipwall}};

Vector<set_bc> HydroState::bc_set = {
    &set_scalar_bc,
    &set_x_vel_bc,
    &set_y_vel_bc,
    &set_z_vel_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc
};

std::string HydroState::tag = "hydro";
bool HydroState::registered = GetStateFactory().Register(HydroState::tag, StateBuilder<HydroState>);

HydroState::HydroState(){}

HydroState::HydroState(const sol::table& def)
{
    num_grow = 0;
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
}

HydroState::~HydroState(){}

#ifdef AMREX_USE_EB
void HydroState::set_eb_bc(const sol::table &bc_def)
{

    std::string bc_type = bc_def.get<std::string>("type");

    if (bc_type == HydroSlipWall::tag) {
        eb_bcs.push_back(std::unique_ptr<HydroBoundaryEB>(new HydroSlipWall(flux_solver.get())));
    } else if (bc_type == HydroNoSlipWall::tag) {
        if (!viscous) {
            Abort("Requested EB bc of type '" + bc_type + "' without defining 'viscosity' for state '" + name + "'");
        }
        eb_bcs.push_back(std::unique_ptr<HydroBoundaryEB>(new HydroNoSlipWall(flux_solver.get(), viscous.get(), bc_def)));
    } else if (bc_type == DirichletWall::tag) {
        eb_bcs.push_back(std::unique_ptr<HydroBoundaryEB>(new DirichletWall(flux_solver.get(), bc_def)));
    } else {
        Abort("Requested EB bc of type '" + bc_type + "' which is not compatible with state '" + name + "'");
    }
}
#endif

void HydroState::set_viscosity()
{

    //
    // viscous terms coefficients
    //

    ClassFactory<HydroViscous> vfact = GetHydroViscousFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string visc = state_def["viscosity"]["type"].get_or<std::string>("");

    viscous = vfact.Build(visc, state_def);

    if (!visc.empty() && !viscous)
        Abort("Invalid viscosity option '"+visc+"'. Options are "+vec2str(vfact.getKeys()));
}

Real HydroState::init_from_number_density(std::map<std::string, Real> data)
{
    Real nd = functions["nd"](data);
    Real alpha = functions["alpha"](data);
    Real m = get_mass(alpha);

    return nd*m;
}

void HydroState::set_udf()
{
    using namespace std::placeholders;

    sol::state& lua = MFP::lua;

    sol::table state_def = lua["states"][name];

    // check if we have 'value' defined
    const sol::table value = state_def["value"].get_or(sol::table());


    if (!value.valid())
        Abort("State "+name+" does not have 'value' defined for initial conditions");

    // get a list of any initialisation functions that need to be called during the run
    sol::table dynamic = state_def["dynamic"].get_or(sol::table());

    bool check, success;

    for (int i = 0; i<prim_names.size(); ++i) {

        std::string comp = prim_names[i];

        // is there a variable with this name?
        success = false;
        check = value[comp].valid();

        // doesn't exist, is there an alternative?
        if (!check) {

            // use number density instead of density
            if (comp.compare("rho") == 0) {

                check = value["nd"].valid();

                if (!check)
                    Abort("State "+name+" does not have 'rho' or 'nd' defined for initial conditions");

                Optional3D1VFunction nd;
                success = get_udf(value["nd"], nd, 0.0);

                functions["nd"] = nd;

                Optional3D1VFunction rho;

                rho.set_func(std::bind(&HydroState::init_from_number_density, this, _1));

                functions[comp] = rho;
            }

        }

        if (!success) {

            Optional3D1VFunction v;

            success = get_udf(value[comp], v, 0.0);

            functions[comp] = v;
        }

        if (dynamic.valid()) {
            for (const auto &d : dynamic) {
                if (d.second.as<std::string>().compare(comp) == 0) {
                    dynamic_functions[i] = &functions[comp];
                }
            }
        }

    }

    return;
}

void HydroState::set_flux()
{

    if (!is_transported()) return;

    ClassFactory<HydroRiemannSolver> rfact = GetHydroRiemannSolverFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string flux = state_def["flux"].get_or<std::string>("null");

    if (flux == "null")
        Abort("Flux option required for state '"+name+"'. Options are "+vec2str(rfact.getKeys()));

    flux_solver = rfact.Build(flux, state_def);

    if (!flux_solver)
        Abort("Invalid flux solver option '"+flux+"'. Options are "+vec2str(rfact.getKeys()));


    return;

}

void HydroState::set_shock_detector()
{

    ClassFactory<HydroShockDetector> sdfact = GetHydroShockDetectorFactory();

    sol::table sd_def = MFP::lua["states"][name]["shock_detector"].get_or(sol::table());

    if (!sd_def.valid()) return;

    sd_def["global_idx"] = global_idx;

    std::string sd_name = sd_def["name"].get_or<std::string>("");

    shock_detector = sdfact.Build(sd_name, sd_def);

    if (!sd_name.empty() && !shock_detector)
        Abort("Invalid shock_detector option '"+sd_name+"'. Options are "+vec2str(sdfact.getKeys()));
}

void HydroState::set_refinement()
{

    ClassFactory<Refinement> rfact = GetHydroRefinementFactory();

    sol::table r_def = MFP::lua["states"][name]["refinement"].get_or(sol::table());

    if (!r_def.valid()) return;

    r_def["global_idx"] = global_idx;

    std::string r_name = r_def["name"].get_or<std::string>("");

    refine = rfact.Build(r_name, r_def);

    if (!r_name.empty() && !refine)
        Abort("Invalid refinement option '"+r_name+"'. Options are "+vec2str(rfact.getKeys()));
}

void HydroState::init_from_lua()
{
    BL_PROFILE("HydroState::init_from_lua");

    EulerianState::init_from_lua();

    sol::state& lua = MFP::lua;

    const sol::table state_def = lua["states"][name];


    //
    // get mass, charge, and density
    //

    expand(state_def["mass"], mass);
    expand(state_def["charge"], charge);
    expand(state_def["gamma"], gamma);

    mass_const = mass[0] == mass[1];
    charge_const = charge[0] == charge[1];
    gamma_const = gamma[0] == gamma[1];

    //
    // user defined functions
    //
    set_udf();

    //
    // viscous terms coefficients
    //
    set_viscosity();


    //
    // domain boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};
    const Vector<std::string>& hydro_var = prim_names;
    const int N = hydro_var.size();

    BoundaryState &bs = boundary_conditions;
    bs.phys_fill_bc.resize(+HydroDef::PrimIdx::NUM);

    for (int ax = 0; ax < AMREX_SPACEDIM; ++ax) {
        for (int lh=0; lh<2; ++lh) {

            std::string side_bc = state_def["bc"][dir_name[ax]][side_name[lh]]["fill_hydro_bc"].get_or<std::string>("outflow");
            int i_side_bc = bc_names.at(side_bc);

            // get any custom values/functions
            for (int j=0; j<N; ++j) {

                if (lh==0) {
                    bs.phys_fill_bc[j].setLo(ax,i_side_bc);
                } else {
                    bs.phys_fill_bc[j].setHi(ax,i_side_bc);
                }

                const sol::object v = state_def["bc"][dir_name[ax]][side_name[lh]][hydro_var[j]].get_or(sol::object());
                Optional3D1VFunction f = get_udf(v);
                bs.set(ax,hydro_var[j],lh,f);

                // special case for inflow condition
                if (i_side_bc == PhysBCType::inflow && !f.is_valid()) {
                    Abort("Setting 'fill_hydro_bc = inflow' requires all primitive variables to be defined, '" + hydro_var[j] + "' is not defined");
                }
            }
#ifdef AMREX_USE_EB
            bool is_symmetry = (i_side_bc == PhysBCType::symmetry) || (i_side_bc == PhysBCType::slipwall) || (i_side_bc == PhysBCType::noslipwall);
            if (lh==0) {
                bs.eb_bc.setLo(ax,is_symmetry ? BCType::reflect_even : BCType::foextrap);
            } else {
                bs.eb_bc.setHi(ax,is_symmetry ? BCType::reflect_even : BCType::foextrap);
            }
#endif
        }
    }

    // check validity of inflow bc
    boundary_conditions.post_init();

    //
    // riemann solver
    //
    set_flux();

    //
    // shock detector
    //
    set_shock_detector();

    //
    // refinement
    //

    set_refinement();

}

void HydroState::variable_setup(Vector<int> periodic)
{

    boundary_conditions.fill_bc.resize(+HydroDef::PrimIdx::NUM);

    for (int icomp=0; icomp < +HydroDef::ConsIdx::NUM; ++icomp) {
        set_bc s = bc_set[icomp]; // the function that sets the bc

        // make sure our periodicity isn't being overwritten
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            if (periodic[d]) {
                boundary_conditions.phys_fill_bc[icomp].setLo(d, PhysBCType::interior);
                boundary_conditions.phys_fill_bc[icomp].setHi(d, PhysBCType::interior);
            }
        }

        // grab the per component BCRec and apply the set bc function to it
        (*s)(boundary_conditions.fill_bc[icomp], boundary_conditions.phys_fill_bc[icomp]);
    }

    Vector<std::string> comp_names(+HydroDef::ConsIdx::NUM);
    for (int icomp=0; icomp < +HydroDef::ConsIdx::NUM; ++icomp) {
        comp_names[icomp] = cons_names[icomp] + "-" + name;
    }

    int ng = num_grow;

#ifdef AMREX_USE_EB
    Interpolater* interp = &eb_cell_cons_interp;
#else
    Interpolater* interp = &cell_cons_interp;
#endif

    bool state_data_extrap = false;
    bool store_in_checkpoint = true;

    DescriptorList& desc_lst = MFP::get_desc_lst();

    data_idx = desc_lst.size();

    desc_lst.addDescriptor(data_idx, IndexType::TheCellType(),
                           StateDescriptor::Point, ng, +HydroDef::ConsIdx::NUM,
                           interp, state_data_extrap,
                           store_in_checkpoint);

    desc_lst.setComponent(
                data_idx, 0, comp_names, boundary_conditions.fill_bc,
                FillBC());


    if (MFP::verbosity >= 1) {
        Print() << str();
    }
}

void HydroState::init_data(MFP* mfp)
{

    const Real* dx = mfp->Geom().CellSize();
    const Real* prob_lo = mfp->Geom().ProbLo();

    MultiFab& S_new = mfp->get_new_data(data_idx);

    Array<Real, +HydroDef::ConsIdx::NUM> U;
    Array<Real, +HydroDef::PrimIdx::NUM> Q;

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        FArrayBox& S_data = S_new[mfi];

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);
        Array4<Real> const& S4 = S_data.array();

#ifdef AMREX_USE_EB
        Array4<const EBCellFlag> const& flag4 = mfp->get_eb_data(global_idx).flags.array(mfi);
#endif

        Real x, y, z;
        for     (int k = lo.z; k <= hi.z; ++k) {
            z = prob_lo[2] + (k + 0.5)*dx[2];
            for   (int j = lo.y; j <= hi.y; ++j) {
                y = prob_lo[1] + (j + 0.5)*dx[1];
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    x = prob_lo[0] + (i + 0.5)*dx[0];

#ifdef AMREX_USE_EB
                    const EBCellFlag &cflag = flag4(i,j,k);

                    if (cflag.isCovered()) {
                        for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                            S4(i,j,k,n) = 0.0;
                        }
                        continue;
                    }
#endif

                    // grab the primitive variables as defined by the user functions
                    for (int icomp=0; icomp<+HydroDef::PrimIdx::NUM; ++icomp) {
                        const std::string& prim_name = prim_names[icomp];
                        const auto& f = functions[prim_name];

                        Q[icomp] = f(x, y, z);

                    }

                    // convert primitive to conserved
                    prim2cons(Q, U);

                    // copy into array
                    for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                        S4(i,j,k,n) = U[n];
                    }
                }
            }
        }
    }
}

Real HydroState::get_density_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    return U[+HydroDef::ConsIdx::Density];
}
Real HydroState::get_density_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q) const
{
    return Q[+HydroDef::PrimIdx::Density];
}


Real HydroState::get_alpha_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    return U[+HydroDef::ConsIdx::Tracer]/U[+HydroDef::ConsIdx::Density];
}
Real HydroState::get_alpha_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q) const
{
    return Q[+HydroDef::PrimIdx::Alpha];
}


Real HydroState::get_mass(Real alpha) const
{
    BL_PROFILE("HydroState::get_mass");

    if (mass_const) return mass[0];

    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    return (m0*m1)/(m0*alpha + m1*(1-alpha));
}

Real HydroState::get_mass(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    BL_PROFILE("HydroState::get_mass");

    if (mass_const) return mass[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_mass(alpha);
}

Real HydroState::get_charge(Real alpha) const
{
    BL_PROFILE("HydroState::get_charge");

    if (charge_const) return charge[0];

    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    Real q0 = charge[0];
    Real q1 = charge[1];

    return (alpha*m0*q1 + (1-alpha)*m1*q0)/(m0*alpha + m1*(1-alpha));
}

Real HydroState::get_charge(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    BL_PROFILE("HydroState::get_charge");

    if (charge_const) return charge[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_charge(alpha);
}

Real HydroState::get_gamma(Real alpha) const
{
    BL_PROFILE("HydroState::get_gamma");

    if (gamma_const) return gamma[0];

    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    Real g0 = gamma[0];
    Real g1 = gamma[1];

    Real cp0 = g0/(m0*(g0-1));
    Real cp1 = g1/(m1*(g1-1));

    Real cv0 = 1/(m0*(g0-1));
    Real cv1 = 1/(m1*(g1-1));

    return ((1-alpha)*cp0 + alpha*cp1)/((1-alpha)*cv0 + alpha*cv1);
}

Real HydroState::get_gamma(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    BL_PROFILE("HydroState::get_gamma_U");

    if (gamma_const) return gamma[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_gamma(alpha);
}

Real HydroState::get_cp(Real alpha) const
{
    BL_PROFILE("HydroState::get_cp");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));


    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    Real g0 = gamma[0];
    Real g1 = gamma[1];

    Real cp0 = g0/(m0*(g0-1));
    Real cp1 = g1/(m1*(g1-1));

    return (1-alpha)*cp0 + alpha*cp1;
}

Real HydroState::get_cp(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    BL_PROFILE("HydroState::get_cp");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));

    Real alpha = get_alpha_from_cons(U);
    return get_cp(alpha);
}


// in place conversion from conserved to primitive
bool HydroState::cons2prim(const Array<Real,+HydroDef::ConsIdx::NUM>& U, Array<Real,+HydroDef::PrimIdx::NUM>& Q) const
{
    BL_PROFILE("HydroState::cons2prim");

    Real rho = U[+HydroDef::ConsIdx::Density];
    Real mx = U[+HydroDef::ConsIdx::Xmom];
    Real my = U[+HydroDef::ConsIdx::Ymom];
    Real mz = U[+HydroDef::ConsIdx::Zmom];
    Real ed = U[+HydroDef::ConsIdx::Eden];
    Real tr = U[+HydroDef::ConsIdx::Tracer];

    Real rhoinv = 1/rho;
    Real u = mx*rhoinv;
    Real v = my*rhoinv;
    Real w = mz*rhoinv;
    Real ke = 0.5*rho*(u*u + v*v + w*w);
    Real alpha = tr*rhoinv;
    Real m = get_mass(alpha);
    Real g = get_gamma(alpha);
    Real p = (ed - ke)*(g - 1);
    Real T = p*rhoinv*m;

    Q[+HydroDef::PrimIdx::Density] = rho;
    Q[+HydroDef::PrimIdx::Xvel] = u;
    Q[+HydroDef::PrimIdx::Yvel] = v;
    Q[+HydroDef::PrimIdx::Zvel] = w;
    Q[+HydroDef::PrimIdx::Prs] = p;
    Q[+HydroDef::PrimIdx::Alpha] = alpha;
    Q[+HydroDef::PrimIdx::Temp] = T;

    return prim_valid(Q);
}

// in-place conversion from primitive to conserved variables
void HydroState::prim2cons(const Array<Real, +HydroDef::PrimIdx::NUM> &Q, Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    BL_PROFILE("HydroState::prim2cons");

    Real rho = Q[+HydroDef::PrimIdx::Density];
    Real u = Q[+HydroDef::PrimIdx::Xvel];
    Real v = Q[+HydroDef::PrimIdx::Yvel];
    Real w = Q[+HydroDef::PrimIdx::Zvel];
    Real p = Q[+HydroDef::PrimIdx::Prs];
    Real alpha = Q[+HydroDef::PrimIdx::Alpha];

    Real mx = u*rho;
    Real my = v*rho;
    Real mz = w*rho;
    Real ke = 0.5*rho*(u*u + v*v + w*w);
    Real tr = alpha*rho;
    Real g = get_gamma(alpha);
    Real ed = p/(g - 1) + ke;

    U[+HydroDef::ConsIdx::Density] = rho;
    U[+HydroDef::ConsIdx::Xmom] = mx;
    U[+HydroDef::ConsIdx::Ymom] = my;
    U[+HydroDef::ConsIdx::Zmom] = mz;
    U[+HydroDef::ConsIdx::Eden] = ed;
    U[+HydroDef::ConsIdx::Tracer] = tr;
    U[+HydroDef::ConsIdx::Shock] = 0.0;

}

bool HydroState::prim_valid(Array<Real,+HydroDef::PrimIdx::NUM>& Q) const
{
    if ((Q[+HydroDef::PrimIdx::Density] <= 0.0) ||  (Q[+HydroDef::PrimIdx::Prs] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

bool HydroState::cons_valid(Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    if ((U[+HydroDef::ConsIdx::Density] <= 0.0) ||  (U[+HydroDef::ConsIdx::Eden] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

Real HydroState::get_energy_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    return U[+HydroDef::ConsIdx::Eden];
}

Real HydroState::get_temperature_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    BL_PROFILE("HydroState::get_temperature_from_cons");

    Real rho = U[+HydroDef::ConsIdx::Density];
    Real mx = U[+HydroDef::ConsIdx::Xmom];
    Real my = U[+HydroDef::ConsIdx::Ymom];
    Real mz = U[+HydroDef::ConsIdx::Zmom];
    Real ed = U[+HydroDef::ConsIdx::Eden];
    Real tr = U[+HydroDef::ConsIdx::Tracer];

    Real rhoinv = 1/rho;
    Real ke = 0.5*rhoinv*(mx*mx + my*my + mz*mz);
    Real alpha = tr*rhoinv;
    Real g = get_gamma(alpha);

    Real prs = (ed - ke)*(g - 1);

    Real m = get_mass(alpha);
    return prs*rhoinv*m;

}

Real HydroState::get_temperature_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q) const
{
    return Q[+HydroDef::PrimIdx::Temp];
}

RealArray HydroState::get_speed_from_cons(const Array<Real,+HydroDef::ConsIdx::NUM>& U) const
{
    BL_PROFILE("HydroState::get_speed_from_cons");

    Real rho = U[+HydroDef::ConsIdx::Density];
    Real mx = U[+HydroDef::ConsIdx::Xmom];
    Real my = U[+HydroDef::ConsIdx::Ymom];
    Real mz = U[+HydroDef::ConsIdx::Zmom];
    Real ed = U[+HydroDef::ConsIdx::Eden];
    Real tr = U[+HydroDef::ConsIdx::Tracer];

    Real rhoinv = 1/rho;

    Real ux = mx*rhoinv;
    Real uy = my*rhoinv;
    Real uz = mz*rhoinv;

    Real kineng = 0.5*rho*(ux*ux + uy*uy + uz*uz);
    Real alpha = tr*rhoinv;
    Real g = get_gamma(alpha);
    Real prs = (ed - kineng)*(g - 1);
    Real a = std::sqrt(g*prs*rhoinv);

    RealArray s = {AMREX_D_DECL(a + std::abs(ux), a + std::abs(uy), a + std::abs(uz))};

    return s;

}

RealArray HydroState::get_speed_from_prim(const Array<Real,+HydroDef::PrimIdx::NUM>& Q) const
{
    BL_PROFILE("HydroState::get_speed_from_prim");

    Real g = get_gamma(Q[+HydroDef::PrimIdx::Alpha]);

    Real a = std::sqrt(g*Q[+HydroDef::PrimIdx::Prs]/Q[+HydroDef::PrimIdx::Density]);

    RealArray s = {AMREX_D_DECL(a + std::abs(Q[+HydroDef::PrimIdx::Xvel]),
                   a + std::abs(Q[+HydroDef::PrimIdx::Yvel]),
                   a + std::abs(Q[+HydroDef::PrimIdx::Zvel]))};


    return s;

}

Vector<std::string> HydroState::get_plot_output_names() const
{
    Vector<std::string> out;
    out.insert(out.end(), cons_names.begin(), cons_names.end());
    out.insert(out.end(), prim_names.begin(), prim_names.end());
    out.push_back("charge");
    out.push_back("mass");
    out.push_back("gamma");
#ifdef AMREX_USE_EB
    out.push_back("vfrac");
#endif
    return out;
}

void HydroState::get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                                 #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                                 #endif
                                 ) const
{
    BL_PROFILE("HydroState::get_state_values");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
    Array4<const Real> const& vf4 = vfrac.array();
#endif

    updated.resize(0);

    // check conserved variables
    std::map<std::string,int> cons_tags;
    for (int i=0; i<+HydroDef::ConsIdx::NUM; ++i) {
        const std::string s = cons_names[i];
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // check primitive variables
    std::map<std::string,int> prim_tags;
    for (int i=0; i<+HydroDef::PrimIdx::NUM; ++i) {
        const std::string s = prim_names[i];
        if (s == cons_names[0]) continue;
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        prim_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // additional variables

    Vector<std::string> other;

    const std::string charge_name = "charge-"+name;
    bool load_charge = out.find(charge_name) != out.end();
    if (load_charge) other.push_back(charge_name);

    const std::string mass_name = "mass-"+name;
    bool load_mass = out.find(mass_name) != out.end();
    if (load_mass) other.push_back(mass_name);

    const std::string gamma_name = "gamma-"+name;
    bool load_gamma = out.find(gamma_name) != out.end();
    if (load_gamma) other.push_back(gamma_name);

#ifdef AMREX_USE_EB
    const std::string vfrac_name = "vfrac-"+name;
    bool load_vfrac = out.find(vfrac_name) != out.end();
    if (load_vfrac) other.push_back(vfrac_name);
#endif

    updated.insert(updated.end(), other.begin(), other.end());

    std::map<std::string,Array4<Real>> out4;
    for (const std::string& s : updated) {
        out[s].resize(box, 1);
        out[s].setVal(0.0);
        out4[s] = out[s].array();
    }

    // temporary storage for retrieving the state data
    Array<Real, +HydroDef::ConsIdx::NUM> S;
    Array<Real, +HydroDef::PrimIdx::NUM> Q;

    Array4<const Real> const& src4 = src.array();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vf4(i,j,k) == 0.0) {
                    continue;
                }
#endif

                for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                    S[n] = src4(i,j,k,n);

                    if (std::isnan(S[n])) {
                        Abort();
                    }

                }

                if (load_charge) out4[charge_name](i,j,k) = get_charge(S);
                if (load_mass)   out4[mass_name](i,j,k)   = get_mass(S);
                if (load_gamma)  out4[gamma_name](i,j,k)  = get_gamma(S);
#ifdef AMREX_USE_EB
                if (load_vfrac)  out4[vfrac_name](i,j,k)  = vf4(i,j,k);
#endif

                if (!cons_tags.empty()) {
                    for (const auto& var : cons_tags) {
                        out4[var.first](i,j,k) = S[var.second];
                    }
                }

                if (!prim_tags.empty()) {
                    cons2prim(S, Q);

                    for (const auto& var : prim_tags) {
                        out4[var.first](i,j,k) = Q[var.second];
                    }
                }
            }
        }
    }


    return;
}

Real HydroState::get_allowed_time_step(MFP* mfp) const
{
    const Real* dx = mfp->Geom().CellSize();

    MultiFab& data = mfp->get_new_data(data_idx);

    Array<Real,+HydroDef::ConsIdx::NUM> U;

    Real max_speed = std::numeric_limits<Real>::min();

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls
        const FArrayBox& vfrac = mfp->get_eb_data(global_idx).volfrac[mfi];

        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();
#endif
        Array4<const Real> const& data4 = data[mfi].array();

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif

                    for (int n = 0; n<+HydroDef::ConsIdx::NUM; ++n) {
                        U[n] = data4(i,j,k,n);
                    }

                    const RealArray speed = get_speed_from_cons(U);


                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        max_speed = std::max(max_speed, speed[d]);
                    }
                }
            }
        }
    }

    Real dt = dx[0]/max_speed;
    for (int i=1; i<AMREX_SPACEDIM; ++i) {
        dt = std::min(dt, dx[i]/max_speed);
    }

    // check for any viscous time step limitation
    if (viscous) {
        dt = std::min(dt, viscous->get_min_dt(mfp));
    }

    return dt;
}

void HydroState::calc_primitives(const Box& box,
                                 FArrayBox& cons,
                                 FArrayBox& prim,
                                 const Real* dx,
                                 const Real t,
                                 const Real* prob_lo
                                 #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                                 #endif
                                 ) const
{
    BL_PROFILE("HydroState::calc_primitives");

    Array<Real, +HydroDef::ConsIdx::NUM> U;
    Array<Real, +HydroDef::PrimIdx::NUM> Q;

    prim.resize(box, +HydroDef::PrimIdx::NUM);

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);
    Array4<Real> const& s4 = cons.array();
    Array4<Real> const& p4 = prim.array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();

    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);
#endif

    Real x, y, z;

    for     (int k = lo.z; k <= hi.z; ++k) {
        z = prob_lo[2] + (k + 0.5)*dx[2];
        for   (int j = lo.y; j <= hi.y; ++j) {
            y = prob_lo[1] + (j + 0.5)*dx[1];
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                x = prob_lo[0] + (i + 0.5)*dx[0];

#ifdef AMREX_USE_EB
                if (vfrac4(i,j,k) == 0.0) {

                    // iterate over all neighbouring cells checking if it has valid data
                    // from these calculate the volume weighted average to populate the
                    // covered cell
                    Real vtot = 0.0;

                    std::fill(U.begin(), U.end(), 0.0);
                    for (const auto& index : grab) {

                        const int ii = i+index[0];
                        const int jj = j+index[1];
                        const int kk = k+index[2];

                        // make sure our stencil is within bounds
                        if ((lo.x > ii) || (ii > hi.x) ||
                                (lo.y > jj) || (jj > hi.y) ||
                                (lo.z > kk) || (kk > hi.z)) continue;

                        const Real vf = vfrac4(ii,jj,kk);
                        if (vf > 0.0) {
                            for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                                U[n] += vf*s4(ii,jj,kk,n);
                            }

                            vtot += vf;
                        }
                    }

                    // if we were close enough to a boundary to have valid data we
                    // average out the volume fraction weighted contributions, otherwise,
                    // fill in the primitives with zeros
                    if (vtot > 0.0) {
                        for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                            U[n] /= vtot;
                        }
                    } else {
                        for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
                            p4(i,j,k,n) = 0.0;
                        }
                        continue;
                    }
                } else {
#endif
                    // grab the conserved variables
                    for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                        U[n] = s4(i,j,k,n);
                    }
#ifdef AMREX_USE_EB
                }
#endif


                // convert to primitive
                cons2prim(U, Q);

                // modify the primitives vector if needed and upload back to
                // the conserved vector
                if (!dynamic_functions.empty()) {
                    for (const auto &f : dynamic_functions) {
                        Q[f.first] = (*f.second)(x, y, z, t);
                    }

                    // copy into primitive
                    for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
                        p4(i,j,k,n) = Q[n];
                    }

                    // convert primitive to conserved
                    prim2cons(Q, U);

                    // copy back into conserved array
                    for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                        s4(i,j,k,n) = U[n];
                    }
                }

                // copy into primitive
                for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
                    p4(i,j,k,n) = Q[n];
                }
            }
        }
    }

    return;
}

void HydroState::update_boundary_cells(const Box& box,
                                       const Geometry &geom,
                                       FArrayBox &prim,
                                       #ifdef AMREX_USE_EB
                                       const FArrayBox& vfrac,
                                       #endif
                                       const Real time) const
{
    BL_PROFILE("HydroState::update_boundary_cells");
    Vector<BoundaryInfo> limits = get_bc_limits(box, geom);

    if (limits.empty())
        return;

    const Box& domain = geom.Domain();
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    Real x, y, z;
    std::map<std::string, Real> Q{{"t", time}};

    Array<int,3> grab;

    Array4<Real> const& p4 = prim.array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();
#endif

    const BCRec &bc = boundary_conditions.where_is_inflow;
    for (const auto &L : limits) {
        if (((L.lo_hi == 0) && (bc.lo(L.dir) == BCType::ext_dir))||
                ((L.lo_hi == 1) && (bc.hi(L.dir) == BCType::ext_dir))) {

            if (L.lo_hi == 0) {
                grab[L.dir] = domain.smallEnd(L.dir);
            } else {
                grab[L.dir] = domain.bigEnd(L.dir);
            }

            for (int k=L.kmin; k<=L.kmax; ++k) {
                z = prob_lo[2] + (k + 0.5)*dx[2];
                Q["z"] = z;
                if (L.dir != 2) {
                    grab[2] = k;
                }
                for (int j=L.jmin; j<=L.jmax; ++j) {
                    y = prob_lo[1] + (j + 0.5)*dx[1];
                    Q["y"] = y;
                    if (L.dir != 1) {
                        grab[1] = j;
                    }
                    for (int i=L.imin; i<=L.imax; ++i) {
                        x = prob_lo[0] + (i + 0.5)*dx[0];
                        Q["x"] = x;
                        if (L.dir != 0) {
                            grab[0] = i;
                        }
#ifdef AMREX_USE_EB
                        if (vfrac4(grab[0],grab[1],grab[2]) == 0.0) {
                            continue;
                        }
#endif

                        // get data from closest internal cell
                        for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
                            Q[prim_names[n]] = p4(grab[0],grab[1],grab[2],n);
                        }

                        // update the primitives from our UDFs, but only those that are valid functions
                        for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
                            const Optional3D1VFunction &f = boundary_conditions.get(L.lo_hi, L.dir, prim_names[n]);
                            if (f.is_valid()) {
                                p4(i,j,k,n) = f(Q);
                            }
                        }
                    }
                }
            }
        }
    }
}

void HydroState::calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                     #ifdef AMREX_USE_EB
                                     ,const EBCellFlagFab &flag
                                     ,const FArrayBox &vfrac
                                     #endif
                                     ) const
{
    BL_PROFILE("HydroState::calc_reconstruction");
    // if we don't want to apply extra limiting on the slopes (forced to 2nd order)
    // we can use the default reconstruction scheme

    // convert pressure
    const Box &pbox = prim.box();
    const Dim3 p_lo = amrex::lbound(pbox);
    const Dim3 p_hi = amrex::ubound(pbox);

    FArrayBox gamma_minus_one(pbox);
    Array4<Real> const& src4 = prim.array();
    Array4<Real> const& gam4 = gamma_minus_one.array();

#ifdef AMREX_USE_EB
    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);

    Array4<const EBCellFlag> const& f4 = flag.array();
    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Vector<Real> stencil(reconstructor->stencil_length);
    int offset = reconstructor->stencil_length/2;
    Array<int,3> stencil_index;
    Array<Real, +HydroDef::PrimIdx::NUM> cell_value, cell_slope;

    Real rho_lo, rho_hi;
    Real alpha_lo, alpha_hi;
    Real abs_phi, phi_scale, coeff_eps;
    Real gam_lo, gam_hi;

    // make sure our arrays for putting lo and hi reconstructed values into
    // are the corect size
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        rlo[d].resize(box, +HydroDef::PrimIdx::NUM);
        rhi[d].resize(box, +HydroDef::PrimIdx::NUM);

#ifdef AMREX_USE_EB
        if (check_eb) {
            rlo[d].copy(prim,box);
            rhi[d].copy(prim,box);
        }
#endif
    }

    // change pressure to internal energy
    for     (int k = p_lo.z; k <= p_hi.z; ++k) {
        for   (int j = p_lo.y; j <= p_hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = p_lo.x; i <= p_hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    continue;
                }
#endif

                gam4(i,j,k) = get_gamma(src4(i,j,k,+HydroDef::PrimIdx::Alpha)) - 1.0;

                src4(i,j,k,+HydroDef::PrimIdx::Prs) /= gam4(i,j,k);

            }
        }
    }

    // now do reconstruction

    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        Array4<Real> const& lo4 = rlo[d].array();
        Array4<Real> const& hi4 = rhi[d].array();

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (check_eb) {

                        // covered cell doesn't need calculating
                        if (f4(i,j,k).isCovered()) {
                            continue;
                        }

                        // cell that references a covered cell doesn't need calculating
                        bool skip = false;
                        stencil_index.fill(0);
                        for (int s=0; s<reconstructor->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            // check if any of the stencil values are from a covered cell
                            if (f4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2]).isCovered()) {
                                skip = true;
                                break;
                            }
                        }

                        if (skip) {
                            continue;
                        }

                    }
#endif

                    // cycle over all components
                    for (int n = 0; n<+HydroDef::PrimIdx::NUM; ++n) {

                        // fill in the stencil along dimension index
                        stencil_index.fill(0);
                        for (int s=0; s<reconstructor->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], n);
                        }

                        // perform reconstruction
                        cell_slope[n] = reconstructor->get_slope(stencil);
                        cell_value[n] = stencil[offset];

                    }


                    // apply corrections to slopes
                    // J. Sci. Comput. (2014) 60:584-611
                    // Robust Finite Volume Schemes for Two-Fluid Plasma Equations

                    Real &rho     = cell_value[+HydroDef::PrimIdx::Density];
                    Real &phi_rho = cell_slope[+HydroDef::PrimIdx::Density];

                    Real &u     = cell_value[+HydroDef::PrimIdx::Xvel];
                    Real &phi_u = cell_slope[+HydroDef::PrimIdx::Xvel];

                    Real &v     = cell_value[+HydroDef::PrimIdx::Yvel];
                    Real &phi_v = cell_slope[+HydroDef::PrimIdx::Yvel];

                    Real &w     = cell_value[+HydroDef::PrimIdx::Zvel];
                    Real &phi_w = cell_slope[+HydroDef::PrimIdx::Zvel];

                    Real &eps     = cell_value[+HydroDef::PrimIdx::Prs];
                    Real &phi_eps = cell_slope[+HydroDef::PrimIdx::Prs];

                    Real &alpha     = cell_value[+HydroDef::PrimIdx::Alpha];
                    Real &phi_alpha = cell_slope[+HydroDef::PrimIdx::Alpha];



                    // correct density slope
                    if (std::abs(phi_rho) > 2*rho) {
                        phi_rho = 2*sign(phi_rho, 0.0)*rho;
                    }

                    // get some face values
                    rho_lo = rho - 0.5*phi_rho;
                    rho_hi = rho + 0.5*phi_rho;

                    alpha_lo = alpha - 0.5*phi_alpha;
                    alpha_hi = alpha + 0.5*phi_alpha;

                    gam_lo = get_gamma(alpha_lo);
                    gam_hi = get_gamma(alpha_hi);

                    abs_phi = phi_u*phi_u + phi_v*phi_v + phi_w*phi_w;

                    // correct velocity slope
                    Real eps_face = eps - 0.5*std::abs(phi_eps);

                    if (eps_face <= 0.0) {
                        // if the reconstructed face value goes non-physical
                        // just set back to first order with zero slope
                        phi_u = 0.0;
                        phi_v = 0.0;
                        phi_w = 0.0;
                        phi_eps = 0.0;
                    } else {
                        coeff_eps = (rho/(rho_lo*rho_hi))*eps_face;
                        if ((0.125*abs_phi) > coeff_eps) {
                            phi_scale = sqrt(abs_phi);
                            coeff_eps = sqrt(8*coeff_eps);
                            phi_u = (phi_u/phi_scale)*coeff_eps;
                            phi_v = (phi_v/phi_scale)*coeff_eps;
                            phi_w = (phi_w/phi_scale)*coeff_eps;
                        }
                        // update eps
                        abs_phi = phi_u*phi_u + phi_v*phi_v + phi_w*phi_w;
                        eps -= (rho_lo*rho_hi/rho)*0.125*abs_phi;
                    }



                    // density
                    lo4(i,j,k,+HydroDef::PrimIdx::Density) = rho_lo;
                    hi4(i,j,k,+HydroDef::PrimIdx::Density) = rho_hi;

                    // x - velocity
                    lo4(i,j,k,+HydroDef::PrimIdx::Xvel) = u - 0.5*(rho_hi/rho)*phi_u;
                    hi4(i,j,k,+HydroDef::PrimIdx::Xvel) = u + 0.5*(rho_lo/rho)*phi_u;

                    // y - velocity
                    lo4(i,j,k,+HydroDef::PrimIdx::Yvel) = v - 0.5*(rho_hi/rho)*phi_v;
                    hi4(i,j,k,+HydroDef::PrimIdx::Yvel) = v + 0.5*(rho_lo/rho)*phi_v;

                    // z - velocity
                    lo4(i,j,k,+HydroDef::PrimIdx::Zvel) = w - 0.5*(rho_hi/rho)*phi_w;
                    hi4(i,j,k,+HydroDef::PrimIdx::Zvel) = w + 0.5*(rho_lo/rho)*phi_w;

                    // epsilon -> pressure
                    lo4(i,j,k,+HydroDef::PrimIdx::Prs) = (eps - 0.5*phi_eps)*(gam_lo - 1.0);
                    hi4(i,j,k,+HydroDef::PrimIdx::Prs) = (eps + 0.5*phi_eps)*(gam_hi - 1.0);

                    Real prs = lo4(i,j,k,+HydroDef::PrimIdx::Prs);

                    // tracer
                    lo4(i,j,k,+HydroDef::PrimIdx::Alpha) = alpha_lo;
                    hi4(i,j,k,+HydroDef::PrimIdx::Alpha) = alpha_hi;

                    // Temperature (calculate from pressure and density)
                    lo4(i,j,k,+HydroDef::PrimIdx::Temp) = lo4(i,j,k,+HydroDef::PrimIdx::Prs)/(rho_lo/get_mass(alpha_lo));
                    hi4(i,j,k,+HydroDef::PrimIdx::Temp) = hi4(i,j,k,+HydroDef::PrimIdx::Prs)/(rho_hi/get_mass(alpha_hi));

                }
            }
        }
    }


    // convert back to pressure
    for     (int k = p_lo.z; k <= p_hi.z; ++k) {
        for   (int j = p_lo.y; j <= p_hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = p_lo.x; i <= p_hi.x; ++i) {
#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif
                src4(i,j,k,+HydroDef::PrimIdx::Prs) *= gam4(i,j,k);



            }
        }
    }

    return;
}



void HydroState::calc_time_averaged_faces(const Box& box,
                                          const FArrayBox &prim,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                          #ifdef AMREX_USE_EB
                                          const EBCellFlagFab& flag,
                                          #endif
                                          const Real* dx,
                                          Real dt) const
{
    BL_PROFILE("HydroState::calc_time_averaged_faces");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& p4 = prim.array();

    Array<Real, +HydroDef::PrimIdx::NUM> Q;

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Real lo_face, centre, hi_face, a6;
    Real dt_2dx;

    const Real ft = 4.0/3.0;

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    continue;
                }
#endif

                for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
                    Q[n] = p4(i,j,k,n);
                }

                const RealArray local_c = get_speed_from_prim(Q);


                for (int d=0; d<AMREX_SPACEDIM; ++d) {

                    Array4<Real> const& lo4 = rlo[d].array();
                    Array4<Real> const& hi4 = rhi[d].array();

                    dt_2dx = 0.5*dt/dx[d];

                    // cycle over all components
                    for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {

                        lo_face = lo4(i,j,k,n);
                        hi_face = hi4(i,j,k,n);
                        centre = p4(i,j,k,n);

                        a6 = 6.0*centre - 3*(lo_face + hi_face);

                        lo4(i,j,k,n) += local_c[d]*dt_2dx*(hi_face - lo_face + (1 - local_c[d]*ft*dt_2dx)*a6);
                        hi4(i,j,k,n) -= local_c[d]*dt_2dx*(hi_face - lo_face - (1 - local_c[d]*ft*dt_2dx)*a6);

                    }
                }
            }
        }
    }

    return;
}

/*
* the data passed to this function is indexed by cell
* but resides at the interface
*/

void HydroState::face_bc(const int dir,
                         Box const& box,
                         const FArrayBox& src,
                         FArrayBox& dest,
                         const Geometry &geom,
                         #ifdef AMREX_USE_EB
                         const EBCellFlagFab& flag,
                         #endif
                         const Real time,
                         const bool do_all) const
{
    BL_PROFILE("HydroState::face_bc");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    const Box& domain = geom.Domain();
    const Dim3 domlo = amrex::lbound(domain);
    const Dim3 domhi = amrex::ubound(domain);

    const Real* prob_lo = geom.ProbLo();
    const Real* dx = geom.CellSize();

    Real x, y, z;
    std::map<std::string, Real> Q{{"t", time}};

    // define the limits of our operations
    Vector<BoundaryInfo> limits;

    BoundaryInfo bi;
    bi.imin = lo.x;
    bi.imax = hi.x;
    bi.jmin = lo.y;
    bi.jmax = hi.y;
    bi.kmin = lo.z;
    bi.kmax = hi.z;



    if (dir == 0) {

        int ilo = domlo.x;
        int ihi = domhi.x + 1; // account for face indexing

        if (lo.x == ilo) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.imin = lo.x;
            b.imax = lo.x;
            b.lo_hi = 0;
        }

        if (hi.x == ihi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.imin = hi.x;
            b.imax = hi.x;
            b.lo_hi = 1;
        }
    }

#if AMREX_SPACEDIM >= 2
    if (dir == 1) {
        int jlo = domlo.y;
        int jhi = domhi.y + 1; // account for face indexing

        if (lo.y == jlo) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.jmin = lo.y;
            b.jmax = lo.y;
            b.lo_hi = 0;
        }

        if (hi.y == jhi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.jmin = hi.y;
            b.jmax = hi.y;
            b.lo_hi = 1;
        }
    }
#endif

#if AMREX_SPACEDIM == 3
    if (dir == 2) {
        int klo = domlo.z;
        int khi = domhi.z + 1; // account for face indexing

        if (lo.z == klo) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.kmin = lo.z;
            b.kmax = lo.z;
            b.lo_hi = 0;
        }

        if (hi.z == khi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.kmin = hi.z;
            b.kmax = hi.z;
            b.lo_hi = 1;
        }
    }
#endif

    Array4<Real> const& d4 = dest.array();
    Array4<const Real> const& s4 = src.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    if (do_all) {
        // first do the usual boundary conditions
        for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
            const BCRec &bc = boundary_conditions.fill_bc[n];

            for (const auto &L : limits) {

                if (
                        ((L.lo_hi == 0) && (bc.lo(dir) == BCType::foextrap))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::foextrap))||

                        ((L.lo_hi == 0) && (bc.lo(dir) == BCType::hoextrap))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::hoextrap))||

                        ((L.lo_hi == 0) && (bc.lo(dir) == BCType::reflect_even))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::reflect_even))
                        ) {
                    for (int k=L.kmin; k<=L.kmax; ++k) {
                        for (int j=L.jmin; j<=L.jmax; ++j) {
                            for (int i=L.imin; i<=L.imax; ++i) {

#ifdef AMREX_USE_EB
                                if (f4(i,j,k).isCovered()) {
                                    continue;
                                }
#endif


                                d4(i,j,k,n) = s4(i,j,k,n);
                            }
                        }
                    }
                }
                if (((L.lo_hi == 0) && (bc.lo(dir) == BCType::reflect_odd))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::reflect_odd))) {
                    for (int k=L.kmin; k<=L.kmax; ++k) {
                        for (int j=L.jmin; j<=L.jmax; ++j) {
                            for (int i=L.imin; i<=L.imax; ++i) {

#ifdef AMREX_USE_EB
                                if (f4(i,j,k).isCovered()) {
                                    continue;
                                }
#endif


                                d4(i,j,k,n) = -s4(i,j,k,n);
                            }
                        }
                    }
                }
            }
        }
    }

    // now do any dirichlet boundaries
    Array<Real,3> offset = {0.5, 0.5, 0.5};
    offset[dir] = 0.0;
    const BCRec &bc = boundary_conditions.where_is_inflow;
    for (const auto &L : limits) {
        if (((L.lo_hi == 0) && (bc.lo(dir) == BCType::ext_dir))||
                ((L.lo_hi == 1) && (bc.hi(dir) == BCType::ext_dir))) {

            for (int k=L.kmin; k<=L.kmax; ++k) {
                z = prob_lo[2] + (k+offset[2])*dx[2];
                Q["z"] = z;
                for (int j=L.jmin; j<=L.jmax; ++j) {
                    y = prob_lo[1] + (j+offset[1])*dx[1];
                    Q["y"] = y;
                    for (int i=L.imin; i<=L.imax; ++i) {
                        x = prob_lo[0] + (i+offset[0])*dx[0];
                        Q["x"] = x;

#ifdef AMREX_USE_EB
                        if (f4(i,j,k).isCovered()) {
                            continue;
                        }
#endif

                        // load data from the other side of the face
                        for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
                            Q[prim_names[n]] = s4(i,j,k,n);
                        }

                        // update the primitives from our UDFs, but only those that are valid functions
                        for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
                            const Optional3D1VFunction &f = boundary_conditions.get(L.lo_hi, dir, prim_names[n]);
                            if (f.is_valid()) {
                                d4(i,j,k,n) = f(Q);
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

// given all of the available face values load the ones expected by the flux calc into a vector
void HydroState::load_state_for_flux(const Array4<const Real> &face,
                                     int i, int j, int k, Array<Real, +HydroDef::FluxIdx::NUM> &S) const
{
    BL_PROFILE("HydroState::load_state_for_flux");


    S[+HydroDef::FluxIdx::Density] = face(i,j,k,+HydroDef::PrimIdx::Density);
    S[+HydroDef::FluxIdx::Xvel] = face(i,j,k,+HydroDef::PrimIdx::Xvel);
    S[+HydroDef::FluxIdx::Yvel] = face(i,j,k,+HydroDef::PrimIdx::Yvel);
    S[+HydroDef::FluxIdx::Zvel] = face(i,j,k,+HydroDef::PrimIdx::Zvel);
    S[+HydroDef::FluxIdx::Prs] = face(i,j,k,+HydroDef::PrimIdx::Prs);
    S[+HydroDef::FluxIdx::Temp] = face(i,j,k,+HydroDef::PrimIdx::Temp);
    S[+HydroDef::FluxIdx::Alpha] = face(i,j,k,+HydroDef::PrimIdx::Alpha);
    S[+HydroDef::FluxIdx::Gamma] = get_gamma(S[+HydroDef::FluxIdx::Alpha]);

    return;
}

void HydroState::calc_fluxes(const Box& box,
                             FArrayBox &cons,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                             Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                             #ifdef AMREX_USE_EB
                             const EBCellFlagFab& flag,
                             #endif
                             const Real *dx,
                             const Real dt) const
{
    BL_PROFILE("HydroState::calc_fluxes");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<int, 3> index;
    Array<Real, +HydroDef::FluxIdx::NUM> L, R;
    Array<Real, +HydroDef::ConsIdx::NUM> F;

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Array4<Real> const& cons4 = cons.array();

    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        FArrayBox& flux = fluxes[d];
        flux.setVal(0.0);
        Array4<Real> const& flux4 = flux.array();

        const Box &fbox = flux.box();
        const Dim3 flo = amrex::lbound(fbox);
        const Dim3 fhi = amrex::ubound(fbox);

        index.fill(0);
        index[d] = 1;

        Array4<const Real> lo4, hi4;

        lo4 = r_lo[d].array();
        hi4 = r_hi[d].array();

        for     (int k = flo.z; k <= fhi.z; ++k) {
            for   (int j = flo.y; j <= fhi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = flo.x; i <= fhi.x; ++i) {

#ifdef AMREX_USE_EB

                    // both cells have to be covered before we skip
                    if (f4(i,j,k).isCovered() && f4(i-index[0],j-index[1],k-index[2]).isCovered())
                        continue;
#endif

                    // get left and right states
                    load_state_for_flux(hi4, i-index[0],j-index[1],k-index[2], L);
                    load_state_for_flux(lo4, i,j,k, R);


                    // rotate the vectors
                    transform_global2local(L, d, flux_vector_idx);
                    transform_global2local(R, d, flux_vector_idx);


                    Real shk = 0.0;

                    if (shock_detector) {
                        shk = shock_detector->solve(L, R);
                        cons4(i-index[0],j-index[1],k-index[2],+HydroDef::ConsIdx::Shock) = std::max(cons4(i-index[0],j-index[1],k-index[2],+HydroDef::ConsIdx::Shock), shk);
                        cons4(i,j,k,+HydroDef::ConsIdx::Shock) = std::max(cons4(i,j,k,+HydroDef::ConsIdx::Shock), shk);
                    }



                    flux_solver->solve(L, R, F, &shk);

                    // rotate the flux back to local frame
                    transform_local2global(F, d, cons_vector_idx);

                    // load the flux into the array
                    for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                        flux4(i,j,k,n) += F[n];
                    }

                }
            }
        }
    }

    return;
}


void HydroState::correct_face_prim(const Box& box,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                                   const Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                   #ifdef AMREX_USE_EB
                                   const EBCellFlagFab &flag,
                                   #endif
                                   const Real *dx,
                                   const Real dt) const
{
    BL_PROFILE("State::correct_face_prim");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<Real, +HydroDef::PrimIdx::NUM> L_prim, R_prim;
    Array<Real, +HydroDef::FluxIdx::NUM> F;
    Array<Real, +HydroDef::ConsIdx::NUM> L_cons, R_cons;

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();

    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    // the below loops selects the alternative dimensions to calculate the corrections from
    // e.g. 1D: d1 = 0, d2 = -
    // e.g. 2D: d1 = 0, d2 = 1
    //          d1 = 1, d2 = 0
    // e.g. 3D: d1 = 0, d2 = 1, 2
    //          d1 = 1, d2 = 2, 0
    //          d1 = 2, d2 = 0, 1

    Array<int, 3> idx1, idx2;

    Real hiF, loF;
    int d1, dd, d2;

    for (d1=0; d1<AMREX_SPACEDIM; ++d1) { // face direction
        idx1.fill(0);
        idx1[d1] = 1;

        Array4<Real> const &lo4 = r_lo[d1].array();
        Array4<Real> const &hi4 = r_hi[d1].array();

        for (dd=1; dd<AMREX_SPACEDIM; ++dd) {
            d2 = (d1+dd)%AMREX_SPACEDIM; // flux direction

            idx2.fill(0);
            idx2[d2] = 1;

            Array4<const Real> const &flux4 = fluxes[d2].array();

            for     (int k = lo.z; k <= hi.z + idx1[2]; ++k) {
                for   (int j = lo.y; j <= hi.y + idx1[1]; ++j) {
                    AMREX_PRAGMA_SIMD
                            for (int i = lo.x; i <= hi.x + idx1[0]; ++i) {

#ifdef AMREX_USE_EB
                        if (check_eb) {
                            // only calculate corrections for faces where the stencil of fluxes is full

                            // a) this interface has a valid flux
                            if (f4(i,j,k).isDisconnected(-idx1[0],-idx1[1],-idx1[2]))
                                continue;

                            // b) right cell has valid fluxes hi and lo
                            if (f4(i,j,k).isDisconnected( idx2[0], idx2[1], idx2[2]))
                                continue;
                            if (f4(i,j,k).isDisconnected(-idx2[0],-idx2[1],-idx2[2]))
                                continue;

                            // c) left cell has valid fluxes hi and lo
                            if (f4(i-idx1[0],j-idx1[1],k-idx1[2]).isDisconnected( idx2[0], idx2[1], idx2[2]))
                                continue;
                            if (f4(i-idx1[0],j-idx1[1],k-idx1[2]).isDisconnected(-idx2[0],-idx2[1],-idx2[2]))
                                continue;
                        }
#endif

                        // get left and right states
                        for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n ) {
                            L_prim[n] = hi4(i-idx1[0],j-idx1[1],k-idx1[2],n);
                            R_prim[n] = lo4(i,j,k,n);
                        }

                        // convert to conserved
                        prim2cons(L_prim, L_cons);
                        prim2cons(R_prim, R_cons);

                        /* example for x-direction in 2D
                                   * L & R = reconstructed x- face values
                                   * hiF & loF = fluxes in y- direction on hi and lo faces
                                   *   relative to the reconstructed values
                                   * diagonal lines represent the contributing factors to L & R
                                   *
                                   * L* = L + 0.5*dt*dy*(loF_L - hiF_L)
                                   * R* = R + 0.5*dt*dy*(loF_R - hiF_R)
                                   *
                                   *
                                   * | - - hiF_L - - | - - hiF_R - - |
                                   * |        \      |     /         |
                                   * |          \    |   /           |
                                   * |             L | R             |
                                   * |          /    |   \           |
                                   * |        /      |     \         |
                                   * | - - loF_L - - | - - loF_R - - |
                                   */

                        // left side
                        for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n ) {
                            loF = flux4(i-idx1[0],j-idx1[1],k-idx1[2],n);
                            hiF = flux4(i-idx1[0]+idx2[0],j-idx1[1]+idx2[1],k-idx1[2]+idx2[2],n);

                            L_cons[n] += 0.5*dt/dx[d2]*(loF - hiF);
                        }

                        // right side
                        for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n ) {
                            loF = flux4(i,j,k,n);
                            hiF = flux4(i+idx2[0],j+idx2[1],k+idx2[2],n);

                            R_cons[n] += 0.5*dt/dx[d2]*(loF - hiF);
                        }

                        // convert to primitive
                        cons2prim(L_cons, L_prim);
                        cons2prim(R_cons, R_prim);

                        // update reconstruction
                        for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n ) {
                            hi4(i-idx1[0],j-idx1[1],k-idx1[2],n) = L_prim[n];
                            lo4(i,j,k,n) = R_prim[n];
                        }
                    }
                }
            }
        }
    }
}

void HydroState::calc_diffusion_terms(const Box& box,
                                      const FArrayBox& prim,
                                      FArrayBox& diff
                                      #ifdef AMREX_USE_EB
                                      ,const EBCellFlagFab& flag
                                      #endif
                                      ) const
{
    BL_PROFILE("HydroState::calc_neutral_diffusion_terms");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& prim4 = prim.array();
    Array4<Real> const& d4 = diff.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Array<Real,+HydroDef::PrimIdx::NUM> Q;
    Real T, mu, kappa;

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                for (int n=0; n<+HydroDef::PrimIdx::NUM; ++n) {
                    Q[n] = prim4(i,j,k,n);
                }

                viscous->get_coeffs(Q, T, mu, kappa);


                d4(i,j,k,+HydroViscous::CoeffIdx::Temp) = T;
                d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) = kappa;
                d4(i,j,k,+HydroViscous::CoeffIdx::Mu) = mu;
            }
        }
    }

    return;
}


void HydroState::calc_viscous_fluxes(const Box& box, Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box& pbox,
                                     const FArrayBox &prim,
                                     #ifdef AMREX_USE_EB
                                     const EBCellFlagFab& flag,
                                     #endif
                                     const Real* dx) const
{
    BL_PROFILE("HydroState::calc_viscous_fluxes");
#ifdef AMREX_USE_EB
    if (flag.getType() != FabType::regular) {
        calc_viscous_fluxes_eb(box, fluxes, pbox, prim, flag, dx);
        return;
    }
#endif

    FArrayBox diff(pbox, +HydroViscous::CoeffIdx::NUM);
    calc_diffusion_terms(pbox,
                         prim,
                         diff
                     #ifdef AMREX_USE_EB
                         ,flag
                     #endif
                         );

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<Real, AMREX_SPACEDIM> dxinv;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dxinv[d] = 1/dx[d];
    }

    Array4<const Real> const& p4 = prim.array();
    Array4<const Real> const& d4 = diff.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Real dudx=0, dudy=0, dudz=0, dvdx=0, dvdy=0, dwdx=0, dwdz=0, divu=0;
    const Real two_thirds = 2/3;

    Real tauxx, tauxy, tauxz, dTdx, muf;

    Array4<Real> const& fluxX = fluxes[0].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x + 1; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdx = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp) - d4(i-1,j,k,+HydroViscous::CoeffIdx::Temp))*dxinv[0];

                dudx = (p4(i,j,k,+HydroDef::PrimIdx::Xvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Xvel))*dxinv[0];
                dvdx = (p4(i,j,k,+HydroDef::PrimIdx::Yvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Yvel))*dxinv[0];
                dwdx = (p4(i,j,k,+HydroDef::PrimIdx::Zvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Zvel))*dxinv[0];

#if AMREX_SPACEDIM >= 2
                dudy = (p4(i,j+1,k,+HydroDef::PrimIdx::Xvel)+p4(i-1,j+1,k,+HydroDef::PrimIdx::Xvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Xvel)-p4(i-1,j-1,k,+HydroDef::PrimIdx::Xvel))*(0.25*dxinv[1]);
                dvdy = (p4(i,j+1,k,+HydroDef::PrimIdx::Yvel)+p4(i-1,j+1,k,+HydroDef::PrimIdx::Yvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Yvel)-p4(i-1,j-1,k,+HydroDef::PrimIdx::Yvel))*(0.25*dxinv[1]);
#endif
#if AMREX_SPACEDIM == 3
                dudz = (p4(i,j,k+1,Xvel)+p4(i-1,j,k+1,Xvel)-p4(i,j,k-1,Xvel)-p4(i-1,j,k-1,Xvel))*(0.25*dxinv[2]);
                dwdz = (p4(i,j,k+1,Zvel)+p4(i-1,j,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i-1,j,k-1,Zvel))*(0.25*dxinv[2]);
#endif
                divu = dudx + dvdy + dwdz;

                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i-1,j,k,+HydroViscous::CoeffIdx::Mu));
                tauxx = muf*(2*dudx-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauxz = muf*(dudz+dwdx);

                fluxX(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxx;
                fluxX(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauxy;
                fluxX(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauxz;
                fluxX(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel) +  p4(i-1,j,k,+HydroDef::PrimIdx::Xvel))*tauxx
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel) + p4(i-1,j,k,+HydroDef::PrimIdx::Yvel))*tauxy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel) + p4(i-1,j,k,+HydroDef::PrimIdx::Zvel))*tauxz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa)+d4(i-1,j,k,+HydroViscous::CoeffIdx::Kappa))*dTdx);

            }
        }
    }

#if AMREX_SPACEDIM >= 2
    Real tauyy, tauyz, dTdy;
    Real dvdz=0, dwdy=0;
    Array4<Real> const& fluxY = fluxes[1].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdy = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp)-d4(i,j-1,k,+HydroViscous::CoeffIdx::Temp))*dxinv[1];
                dudy = (p4(i,j,k,+HydroDef::PrimIdx::Xvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Xvel))*dxinv[1];
                dvdy = (p4(i,j,k,+HydroDef::PrimIdx::Yvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Yvel))*dxinv[1];
                dwdy = (p4(i,j,k,+HydroDef::PrimIdx::Zvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Zvel))*dxinv[1];
                dudx = (p4(i+1,j,k,+HydroDef::PrimIdx::Xvel)+p4(i+1,j-1,k,+HydroDef::PrimIdx::Xvel)-p4(i-1,j,k,+HydroDef::PrimIdx::Xvel)-p4(i-1,j-1,k,+HydroDef::PrimIdx::Xvel))*(0.25*dxinv[0]);
                dvdx = (p4(i+1,j,k,+HydroDef::PrimIdx::Yvel)+p4(i+1,j-1,k,+HydroDef::PrimIdx::Yvel)-p4(i-1,j,k,+HydroDef::PrimIdx::Yvel)-p4(i-1,j-1,k,+HydroDef::PrimIdx::Yvel))*(0.25*dxinv[0]);
#if AMREX_SPACEDIM == 3
                dvdz = (p4(i,j,k+1,Yvel)+p4(i,j-1,k+1,Yvel)-p4(i,j,k-1,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[2]);
                dwdz = (p4(i,j,k+1,Zvel)+p4(i,j-1,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i,j-1,k-1,Zvel))*(0.25*dxinv[2]);
#endif
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i,j-1,k,+HydroViscous::CoeffIdx::Mu));
                tauyy = muf*(2*dvdy-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauyz = muf*(dwdy+dvdz);

                fluxY(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxy;
                fluxY(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauyy;
                fluxY(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauyz;
                fluxY(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Xvel))*tauxy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Yvel))*tauyy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Zvel))*tauyz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) + d4(i,j-1,k,+HydroViscous::CoeffIdx::Kappa))*dTdy);

            }
        }
    }


#endif
#if AMREX_SPACEDIM == 3
    Real tauzz, dTdz;
    Array4<Real> const& fluxZ = fluxes[2].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdz = (d4(i,j,k,iTemp)-d4(i,j,k-1,iTemp))*dxinv[2];
                dudz = (p4(i,j,k,Xvel)-p4(i,j,k-1,Xvel))*dxinv[2];
                dvdz = (p4(i,j,k,Yvel)-p4(i,j,k-1,Yvel))*dxinv[2];
                dwdz = (p4(i,j,k,Zvel)-p4(i,j,k-1,Zvel))*dxinv[2];
                dudx = (p4(i+1,j,k,Xvel)+p4(i+1,j,k-1,Xvel)-p4(i-1,j,k,Xvel)-p4(i-1,j,k-1,Xvel))*(0.25*dxinv[0]);
                dwdx = (p4(i+1,j,k,Zvel)+p4(i+1,j,k-1,Zvel)-p4(i-1,j,k,Zvel)-p4(i-1,j,k-1,Zvel))*(0.25*dxinv[0]);
                dvdy = (p4(i,j+1,k,Yvel)+p4(i,j+1,k-1,Yvel)-p4(i,j-1,k,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[1]);
                dwdy = (p4(i,j+1,k,Zvel)+p4(i,j+1,k-1,Zvel)-p4(i,j-1,k,Zvel)-p4(i,j-1,k-1,Zvel))*(0.25*dxinv[1]);
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,iMu)+d4(i,j,k-1,iMu));
                tauxz = muf*(dudz+dwdx);
                tauyz = muf*(dvdz+dwdy);
                tauzz = muf*(2.*dwdz-two_thirds*divu);

                fluxZ(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauyz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauzz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Xvel))*tauxz
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Yvel))*tauyz
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Zvel))*tauzz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) +d4(i,j,k-1,+HydroViscous::CoeffIdx::Kappa))*dTdz);

            }
        }
    }

#endif


    return;
}

#ifdef AMREX_USE_EB
void HydroState::calc_viscous_fluxes_eb(const Box& box, Array<FArrayBox,
                                        AMREX_SPACEDIM> &fluxes,
                                        const Box& pbox,
                                        const FArrayBox &prim,
                                        #ifdef AMREX_USE_EB
                                        const EBCellFlagFab& flag,
                                        #endif
                                        const Real* dx) const
{
    BL_PROFILE("HydroState::calc_neutral_viscous_fluxes_eb");
    FArrayBox diff(pbox, +HydroViscous::CoeffIdx::NUM);
    calc_diffusion_terms(pbox,
                         prim,
                         diff
                     #ifdef AMREX_USE_EB
                         ,flag
                     #endif
                         );

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<Real, AMREX_SPACEDIM> dxinv;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dxinv[d] = 1/dx[d];
    }

    Array4<const Real> const& p4 = prim.array();
    Array4<const Real> const& d4 = diff.array();

    Array4<const EBCellFlag> const& f4 = flag.array();

    Real dudx=0, dudy=0, dudz=0, dvdx=0, dvdy=0, dwdx=0, dwdz=0, divu=0;
    const Real two_thirds = 2/3;

    const Array<Real,3>  weights = {0.0, 1.0, 0.5};
    Real whi, wlo;

    Real tauxx, tauxy, tauxz, dTdx, muf;


    // X - direction
    Array4<Real> const& fluxX = fluxes[0].array();
    for     (int k = lo.z-AMREX_D_PICK(0,0,1); k <= hi.z+AMREX_D_PICK(0,0,1); ++k) {
        for   (int j = lo.y-AMREX_D_PICK(0,1,1); j <= hi.y+AMREX_D_PICK(0,1,1); ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x + 1; ++i) {

                bool covered = f4(i,j,k).isCovered();
                bool connected = f4(i,j,k).isConnected(-1,0,0);
                bool other_covered = f4(i-1,j,k).isCovered();

                // only calculate fluxes for fluid cells and between cells that are connected
                if (covered || other_covered || !connected)
                    continue;

                dTdx = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp) - d4(i-1,j,k,+HydroViscous::CoeffIdx::Temp))*dxinv[0];

                dudx = (p4(i,j,k,+HydroDef::PrimIdx::Xvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Xvel))*dxinv[0];
                dvdx = (p4(i,j,k,+HydroDef::PrimIdx::Yvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Yvel))*dxinv[0];
                dwdx = (p4(i,j,k,+HydroDef::PrimIdx::Zvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Zvel))*dxinv[0];

#if AMREX_SPACEDIM >= 2

                const int jhip = j + (int)f4(i,j,k).isConnected(0, 1,0);
                const int jhim = j - (int)f4(i,j,k).isConnected(0,-1,0);
                const int jlop = j + (int)f4(i-1,j,k).isConnected(0, 1,0);
                const int jlom = j - (int)f4(i-1,j,k).isConnected(0,-1,0);
                whi = weights[jhip-jhim];
                wlo = weights[jlop-jlom];
                dudy = (0.5*dxinv[1]) * ((p4(i  ,jhip,k,+HydroDef::PrimIdx::Xvel)-p4(i  ,jhim,k,+HydroDef::PrimIdx::Xvel))*whi+(p4(i-1,jlop,k,+HydroDef::PrimIdx::Xvel)-p4(i-1,jlom,k,+HydroDef::PrimIdx::Xvel))*wlo);
                dvdy = (0.50*dxinv[1]) * ((p4(i  ,jhip,k,+HydroDef::PrimIdx::Yvel)-p4(i  ,jhim,k,+HydroDef::PrimIdx::Yvel))*whi+(p4(i-1,jlop,k,+HydroDef::PrimIdx::Yvel)-p4(i-1,jlom,k,+HydroDef::PrimIdx::Yvel))*wlo);

#endif
#if AMREX_SPACEDIM == 3

                const int khip = k + (int)f4(i,j,k).isConnected(0,0, 1);
                const int khim = k - (int)f4(i,j,k).isConnected(0,0,-1);
                const int klop = k + (int)f4(i-1,j,k).isConnected(0,0, 1);
                const int klom = k - (int)f4(i-1,j,k).isConnected(0,0,-1);
                whi = weights[khip-khim];
                wlo = weights[klop-klom];
                dudz = (0.5*dxinv[2]) * ((p4(i  ,j,khip,+HydroDef::PrimIdx::Xvel)-p4(i  ,j,khim,+HydroDef::PrimIdx::Xvel))*whi + (p4(i-1,j,klop,+HydroDef::PrimIdx::Xvel)-p4(i-1,j,klom,+HydroDef::PrimIdx::Xvel))*wlo);
                dwdz = (0.5*dxinv[2]) * ((p4(i  ,j,khip,+HydroDef::PrimIdx::Zvel)-p4(i  ,j,khim,+HydroDef::PrimIdx::Zvel))*whi + (p4(i-1,j,klop,+HydroDef::PrimIdx::Zvel)-p4(i-1,j,klom,+HydroDef::PrimIdx::Zvel))*wlo);

#endif
                divu = dudx + dvdy + dwdz;

                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i-1,j,k,+HydroViscous::CoeffIdx::Mu));
                tauxx = muf*(2*dudx-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauxz = muf*(dudz+dwdx);

                fluxX(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxx;
                fluxX(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauxy;
                fluxX(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauxz;
                fluxX(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel) +  p4(i-1,j,k,+HydroDef::PrimIdx::Xvel))*tauxx
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel) + p4(i-1,j,k,+HydroDef::PrimIdx::Yvel))*tauxy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel) + p4(i-1,j,k,+HydroDef::PrimIdx::Zvel))*tauxz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa)+d4(i-1,j,k,+HydroViscous::CoeffIdx::Kappa))*dTdx);
            }
        }
    }

    // Y - direction
#if AMREX_SPACEDIM >= 2
    Real tauyy, tauyz, dTdy;
    Real dvdz=0, dwdy=0;
    Array4<Real> const& fluxY = fluxes[1].array();
    for     (int k = lo.z-AMREX_D_PICK(0,0,1); k <= hi.z+AMREX_D_PICK(0,0,1); ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-1; i <= hi.x+1; ++i) {

                bool covered = f4(i,j,k).isCovered();
                bool connected = f4(i,j,k).isConnected(0,-1,0);
                bool other_covered = f4(i,j-1,k).isCovered();

                // only calculate fluxes for fluid cells and between cells that are connected
                if (covered || other_covered || !connected)
                    continue;

                dTdy = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp)-d4(i,j-1,k,+HydroViscous::CoeffIdx::Temp))*dxinv[1];
                dudy = (p4(i,j,k,+HydroDef::PrimIdx::Xvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Xvel))*dxinv[1];
                dvdy = (p4(i,j,k,+HydroDef::PrimIdx::Yvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Yvel))*dxinv[1];
                dwdy = (p4(i,j,k,+HydroDef::PrimIdx::Zvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Zvel))*dxinv[1];

                const int ihip = i + (int)f4(i,j,  k).isConnected( 1,0,0);
                const int ihim = i - (int)f4(i,j,  k).isConnected(-1,0,0);
                const int ilop = i + (int)f4(i,j-1,k).isConnected( 1,0,0);
                const int ilom = i - (int)f4(i,j-1,k).isConnected(-1,0,0);
                whi = weights[ihip-ihim];
                wlo = weights[ilop-ilom];

                dudx = (0.5*dxinv[0]) * ((p4(ihip,j  ,k,+HydroDef::PrimIdx::Xvel)-p4(ihim,j  ,k,+HydroDef::PrimIdx::Xvel))*whi + (p4(ilop,j-1,k,+HydroDef::PrimIdx::Xvel)-p4(ilom,j-1,k,+HydroDef::PrimIdx::Xvel))*wlo);
                dvdx = (0.5*dxinv[0]) * ((p4(ihip,j  ,k,+HydroDef::PrimIdx::Yvel)-p4(ihim,j  ,k,+HydroDef::PrimIdx::Yvel))*whi + (p4(ilop,j-1,k,+HydroDef::PrimIdx::Yvel)-p4(ilom,j-1,k,+HydroDef::PrimIdx::Yvel))*wlo);

#if AMREX_SPACEDIM == 3

                const int khip = k + (int)f4(i,j,  k).isConnected(0,0, 1);
                const int khim = k - (int)f4(i,j,  k).isConnected(0,0,-1);
                const int klop = k + (int)f4(i,j-1,k).isConnected(0,0, 1);
                const int klom = k - (int)f4(i,j-1,k).isConnected(0,0,-1);
                whi = weights[khip-khim];
                wlo = weights[klop-klom];

                dvdz = (0.5*dxinv[2]) * ((p4(i,j  ,khip,+HydroDef::PrimIdx::Yvel)-p4(i,j  ,khim,+HydroDef::PrimIdx::Yvel))*whi + (p4(i,j-1,klop,+HydroDef::PrimIdx::Yvel)-p4(i,j-1,klom,+HydroDef::PrimIdx::Yvel))*wlo);
                dwdz = (0.5*dxinv[2]) * ((p4(i,j  ,khip,+HydroDef::PrimIdx::Zvel)-p4(i,j  ,khim,+HydroDef::PrimIdx::Zvel))*whi + (p4(i,j-1,klop,+HydroDef::PrimIdx::Zvel)-p4(i,j-1,klom,+HydroDef::PrimIdx::Zvel))*wlo);

#endif
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i,j-1,k,+HydroViscous::CoeffIdx::Mu));
                tauyy = muf*(2*dvdy-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauyz = muf*(dwdy+dvdz);

                fluxY(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxy;
                fluxY(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauyy;
                fluxY(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauyz;
                fluxY(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Xvel))*tauxy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Yvel))*tauyy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Zvel))*tauyz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) + d4(i,j-1,k,+HydroViscous::CoeffIdx::Kappa))*dTdy);

            }
        }
    }
#endif

    // Z - direction
#if AMREX_SPACEDIM == 3
    Real tauzz, dTdz;
    Array4<Real> const& fluxZ = fluxes[2].array();
    for     (int k = lo.z; k <= hi.z+1; ++k) {
        for   (int j = lo.y-1; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-1; i <= hi.x+1; ++i) {

                bool covered = f4(i,j,k).isCovered();
                bool connected = f4(i,j,k).isConnected(0,0,-1);
                bool other_covered = f4(i,j,k-1).isCovered();

                // only calculate fluxes for fluid cells and between cells that are connected
                if (covered || other_covered || !connected)
                    continue;

                dTdz = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp)-d4(i,j,k-1,+HydroViscous::CoeffIdx::Temp))*dxinv[2];
                dudz = (p4(i,j,k,+HydroDef::PrimIdx::Xvel)-p4(i,j,k-1,+HydroDef::PrimIdx::Xvel))*dxinv[2];
                dvdz = (p4(i,j,k,+HydroDef::PrimIdx::Yvel)-p4(i,j,k-1,+HydroDef::PrimIdx::Yvel))*dxinv[2];
                dwdz = (p4(i,j,k,+HydroDef::PrimIdx::Zvel)-p4(i,j,k-1,+HydroDef::PrimIdx::Zvel))*dxinv[2];

                const int ihip = i + (int)f4(i,j,k  ).isConnected( 1,0,0);
                const int ihim = i - (int)f4(i,j,k  ).isConnected(-1,0,0);
                const int ilop = i + (int)f4(i,j,k-1).isConnected( 1,0,0);
                const int ilom = i - (int)f4(i,j,k-1).isConnected(-1,0,0);
                whi = weights[ihip-ihim];
                wlo = weights[ilop-ilom];

                dudx = (0.5*dxinv[0]) * ((p4(ihip,j,k  ,+HydroDef::PrimIdx::Xvel)-p4(ihim,j,k  ,+HydroDef::PrimIdx::Xvel))*whi + (p4(ilop,j,k-1,+HydroDef::PrimIdx::Xvel)-p4(ilom,j,k-1,+HydroDef::PrimIdx::Xvel))*wlo);
                dwdx = (0.5*dxinv[0]) * ((p4(ihip,j,k  ,+HydroDef::PrimIdx::Zvel)-p4(ihim,j,k  ,+HydroDef::PrimIdx::Zvel))*whi + (p4(ilop,j,k-1,+HydroDef::PrimIdx::Zvel)-p4(ilom,j,k-1,+HydroDef::PrimIdx::Zvel))*wlo);

                const int jhip = j + (int)f4(i,j,k  ).isConnected(0 ,1,0);
                const int jhim = j - (int)f4(i,j,k  ).isConnected(0,-1,0);
                const int jlop = j + (int)f4(i,j,k-1).isConnected(0 ,1,0);
                const int jlom = j - (int)f4(i,j,k-1).isConnected(0,-1,0);
                whi = weights[jhip-jhim];
                wlo = weights[jlop-jlom];

                dvdy = (0.5*dxinv[1]) * ((p4(i,jhip,k  ,+HydroDef::PrimIdx::Yvel)-p4(i,jhim,k  ,+HydroDef::PrimIdx::Yvel))*whi + (p4(i,jlop,k-1,+HydroDef::PrimIdx::Yvel)-p4(i,jlom,k-1,+HydroDef::PrimIdx::Yvel))*wlo);
                dwdy = (0.5*dxinv[1]) * ((p4(i,jhip,k  ,+HydroDef::PrimIdx::Zvel)-p4(i,jhim,k  ,+HydroDef::PrimIdx::Zvel))*whi + (p4(i,jlop,k-1,+HydroDef::PrimIdx::Zvel)-p4(i,jlom,k-1,+HydroDef::PrimIdx::Zvel))*wlo);

                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i,j,k-1,+HydroViscous::CoeffIdx::Mu));
                tauxz = muf*(dudz+dwdx);
                tauyz = muf*(dvdz+dwdy);
                tauzz = muf*(2.*dwdz-two_thirds*divu);

                fluxZ(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauyz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauzz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Xvel))*tauxz
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Yvel))*tauyz
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Zvel))*tauzz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) +d4(i,j,k-1,+HydroViscous::CoeffIdx::Kappa))*dTdz);

            }
        }
    }

#endif

}


void HydroState::calc_wall_fluxes(const Box& box,
                                  const FArrayBox &prim,
                                  Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                  const EBCellFlagFab& flag,
                                  const CutFab &bc_idx,
                                  const FArrayBox& bcent,
                                  const FArrayBox &bnorm,
                                  const Array<const FArrayBox*, AMREX_SPACEDIM> &afrac,
                                  const Real *dx,
                                  const Real dt) const
{
    BL_PROFILE("HydroState::calc_wall_fluxes");

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& p4 = prim.array();
    Array4<const Real> const& bcent4 = bcent.array();
    Array4<const Real> const& bnorm4 = bnorm.array();

    Array<Array4<Real>,AMREX_SPACEDIM> flux4;
    Array<Array4<const Real>,AMREX_SPACEDIM> afrac4;

    Array<Array<Real,+HydroDef::ConsIdx::NUM>,AMREX_SPACEDIM> wall_flux;

    Array4<const EBCellFlag> const& flag4 = flag.array();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
        afrac4[d] = afrac[d]->array();

        // zero out the flux accumulator
        fluxes[d].setVal(0.0);
    }

    const Array4<const Real>& bc_idx4 = bc_idx.array();

    Array<Real,+HydroDef::FluxIdx::NUM> cell_state;

    Array<Array<Real,3>,3> wall_coord = {{{0,0,0},{0,0,0},{0,0,0}}};
    Array<Real,AMREX_SPACEDIM> wall_centre;

    for (int k = lo.z-AMREX_D_PICK(0,0,2); k <= hi.z+AMREX_D_PICK(0,0,2); ++k) {
        for (int j = lo.y-AMREX_D_PICK(0,2,2); j <= hi.y+AMREX_D_PICK(0,2,2); ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-2; i <= hi.x+2; ++i) {

                const EBCellFlag &cflag = flag4(i,j,k);

                if (cflag.isSingleValued()) {


                    // grab a vector of the local state
                    load_state_for_flux(p4, i, j, k, cell_state);

                    for (int d=0; d<AMREX_SPACEDIM; ++d) {

                        // get the wall normal
                        wall_coord[0][d] = bnorm4(i,j,k,d);

                        // get the centre of the wall
                        wall_centre[d] = bcent4(i,j,k,d);
                    }

                    // get a local coordinate system with x- aligned with the wall normal
                    expand_coord(wall_coord);

                    // the boundary condition
                    const int ebi = (int)nearbyint(bc_idx4(i,j,k));
                    const HydroBoundaryEB& bc = *eb_bcs[ebi];

                    // calculate the wall flux
                    bc.solve(wall_coord, wall_centre, cell_state, p4, i, j, k, dx, wall_flux);

                    // load the flux into the fab
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                            flux4[d](i,j,k,n) += wall_flux[d][n];
                        }
                    }
                }
            }
        }
    }
    return;
}

#endif

void HydroState::write_info(nlohmann::json &js) const
{

    EulerianState::write_info(js);

    // write out stuff that is common to all states

    js["type"] = tag;
    js["type_idx"] = +get_type();

    js["mass"] = mass;
    js["charge"] = charge;
    js["gamma"] = gamma;

    if (viscous) {

        auto& grp = js["viscosity"];

        grp["type"] = viscous->get_tag();

        const auto coeffs = viscous->get_refs();

        for (const auto& cf : coeffs) {
            grp[cf.first] = cf.second;
        }
    }

}
