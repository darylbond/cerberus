#include "MFP_hydro.H"
#include "MFP_lua.H"
#include "MFP_hydro_bc.H"
#include "MFP.H"

Vector<std::string> HydroState::cons_names = {
    "rho",
    "x_mom",
    "y_mom",
    "z_mom",
    "nrg",
    "tracer"
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

Vector<set_bc> HydroState::bc_set = {
    &set_scalar_bc,
    &set_x_vel_bc,
    &set_y_vel_bc,
    &set_z_vel_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc
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

std::string HydroState::tag = "hydro";
bool HydroState::registered = GetStateFactory().Register(HydroState::tag, StateBuilder<HydroState>);

HydroState::HydroState() {}

HydroState::HydroState(const sol::table &def)
{
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
}

HydroState::~HydroState(){}


void HydroState::set_viscosity()
{

    //
    // viscous terms coefficients
    //

    ClassFactory<HydroViscous> vfact = GetViscousFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string visc = state_def["viscosity"]["type"].get_or<std::string>("");

    viscous = vfact.Build(visc, state_def);

    if (!visc.empty() && !viscous)
        Abort("Invalid viscosity option '"+visc+"'. Options are "+vec2str(vfact.getKeys()));

}

#ifdef AMREX_USE_EB
void HydroState::set_eb_bc(const sol::table &bc_def)
{

    std::string bc_type = bc_def.get<std::string>("type");

    if (bc_type == HydroSlipWall::tag) {
        eb_bcs.push_back(std::unique_ptr<BoundaryEB>(new HydroSlipWall(flux_solver.get())));
    } else if (bc_type == HydroNoSlipWall::tag) {
        if (!viscous) {
            Abort("Requested EB bc of type '" + bc_type + "' without defining 'viscosity' for state '" + name + "'");
        }
        eb_bcs.push_back(std::unique_ptr<BoundaryEB>(new HydroNoSlipWall(flux_solver.get(), viscous.get(), bc_def)));
    } else if (bc_type == DirichletWall::tag) {
        eb_bcs.push_back(std::unique_ptr<BoundaryEB>(new DirichletWall(flux_solver.get(), prim_names, prim_vector_idx, bc_def)));
    } else {
        Abort("Requested EB bc of type '" + bc_type + "' which is not compatible with state '" + name + "'");
    }
}
#endif


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

    sol::state& lua = MFP::lua;

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
    // shock detector threshold
    //
    set_shock_detector();


    //
    // positivity
    //
    enforce_positivity = state_def["enforce_positivity"].get_or(0);
    extra_slope_limits = state_def["extra_slope_limits"].get_or(1);

}


Real HydroState::get_density_from_cons(const Array<Real,+ConsIdx::NUM>& U)
{
    return U[+ConsIdx::Density];
}
Real HydroState::get_density_from_prim(const Array<Real,+PrimIdx::NUM>& Q)
{
    return Q[+PrimIdx::Density];
}


Real HydroState::get_alpha_from_cons(const Array<Real,+ConsIdx::NUM>& U)
{
    return U[+ConsIdx::Tracer]/U[+ConsIdx::Density];
}
Real HydroState::get_alpha_from_prim(const Array<Real,+PrimIdx::NUM>& Q)
{
    return Q[+PrimIdx::Alpha];
}


Real HydroState::get_mass(Real alpha)
{
    BL_PROFILE("HydroState::get_mass");

    if (mass_const) return mass[0];

    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    return (m0*m1)/(m0*alpha + m1*(1-alpha));
}

Real HydroState::get_mass(const Array<Real,+ConsIdx::NUM>& U)
{
    BL_PROFILE("HydroState::get_mass");

    if (mass_const) return mass[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_mass(alpha);
}

Real HydroState::get_charge(Real alpha)
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

Real HydroState::get_charge(const Array<Real,+ConsIdx::NUM>& U)
{
    BL_PROFILE("HydroState::get_charge");

    if (charge_const) return charge[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_charge(alpha);
}

Real HydroState::get_gamma(Real alpha)
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

Real HydroState::get_gamma(const Array<Real,+ConsIdx::NUM>& U)
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

Real HydroState::get_cp(const Array<Real,+ConsIdx::NUM>& U)
{
    BL_PROFILE("HydroState::get_cp");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));

    Real alpha = get_alpha_from_cons(U);
    return get_cp(alpha);
}


// in place conversion from conserved to primitive
bool HydroState::cons2prim(Vector<Real>& U, Vector<Real>& Q)
{
    BL_PROFILE("HydroState::cons2prim");

    Real rho = U[+ConsIdx::Density];
    Real mx = U[+ConsIdx::Xmom];
    Real my = U[+ConsIdx::Ymom];
    Real mz = U[+ConsIdx::Zmom];
    Real ed = U[+ConsIdx::Eden];
    Real tr = U[+ConsIdx::Tracer];

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

    Q[+PrimIdx::Density] = rho;
    Q[+PrimIdx::Xvel] = u;
    Q[+PrimIdx::Yvel] = v;
    Q[+PrimIdx::Zvel] = w;
    Q[+PrimIdx::Prs] = p;
    Q[+PrimIdx::Alpha] = alpha;
    Q[+PrimIdx::Temp] = T;

    return prim_valid(Q);
}

// in-place conversion from primitive to conserved variables
void HydroState::prim2cons(Vector<Real>& Q, Vector<Real>& U)
{
    BL_PROFILE("HydroState::prim2cons");

    Real rho = Q[+PrimIdx::Density];
    Real u = Q[+PrimIdx::Xvel];
    Real v = Q[+PrimIdx::Yvel];
    Real w = Q[+PrimIdx::Zvel];
    Real p = Q[+PrimIdx::Prs];
    Real alpha = Q[+PrimIdx::Alpha];

    Real mx = u*rho;
    Real my = v*rho;
    Real mz = w*rho;
    Real ke = 0.5*rho*(u*u + v*v + w*w);
    Real tr = alpha*rho;
    Real g = get_gamma(alpha);
    Real ed = p/(g - 1) + ke;

    U[+ConsIdx::Density] = rho;
    U[+ConsIdx::Xmom] = mx;
    U[+ConsIdx::Ymom] = my;
    U[+ConsIdx::Zmom] = mz;
    U[+ConsIdx::Eden] = ed;
    U[+ConsIdx::Tracer] = tr;

}


bool HydroState::prim_valid(Array<Real,+PrimIdx::NUM>& Q)
{
    if ((Q[+PrimIdx::Density] <= 0.0) ||  (Q[+PrimIdx::Prs] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

bool HydroState::cons_valid(Array<Real,+ConsIdx::NUM>& U)
{
    if ((U[+ConsIdx::Density] <= 0.0) ||  (U[+ConsIdx::Eden] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

Real HydroState::get_energy_from_cons(const Array<Real,+ConsIdx::NUM>& U)
{
    return U[+ConsIdx::Eden];
}

Real HydroState::get_temperature_from_cons(const Array<Real,+ConsIdx::NUM>& U)
{
    BL_PROFILE("HydroState::get_temperature_from_cons");

    Real rho = U[+ConsIdx::Density];
    Real mx = U[+ConsIdx::Xmom];
    Real my = U[+ConsIdx::Ymom];
    Real mz = U[+ConsIdx::Zmom];
    Real ed = U[+ConsIdx::Eden];
    Real tr = U[+ConsIdx::Tracer];

    Real rhoinv = 1/rho;
    Real ke = 0.5*rhoinv*(mx*mx + my*my + mz*mz);
    Real alpha = tr*rhoinv;
    Real g = get_gamma(alpha);

    Real prs = (ed - ke)*(g - 1);

    Real m = get_mass(alpha);
    return prs*rhoinv*m;

}

Real HydroState::get_temperature_from_prim(const Array<Real,+PrimIdx::NUM>& Q)
{
    return Q[+PrimIdx::Temp];
}

RealArray HydroState::get_speed_from_cons(const Array<Real,+ConsIdx::NUM>& U)
{
    BL_PROFILE("HydroState::get_speed_from_cons");

    Real rho = U[+ConsIdx::Density];
    Real mx = U[+ConsIdx::Xmom];
    Real my = U[+ConsIdx::Ymom];
    Real mz = U[+ConsIdx::Zmom];
    Real ed = U[+ConsIdx::Eden];
    Real tr = U[+ConsIdx::Tracer];

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

RealArray HydroState::get_speed_from_prim(const Vector<Real>& Q)
{
    BL_PROFILE("HydroState::get_speed_from_prim");

    Real g = get_gamma(Q[+PrimIdx::Alpha]);

    Real a = std::sqrt(g*Q[+PrimIdx::Prs]/Q[+PrimIdx::Density]);

    RealArray s = {AMREX_D_DECL(a + std::abs(Q[+PrimIdx::Xvel]),
                   a + std::abs(Q[+PrimIdx::Yvel]),
                   a + std::abs(Q[+PrimIdx::Zvel]))};


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
    out.push_back("vfrac");
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
    for (int i=0; i<+ConsIdx::NUM; ++i) {
        const std::string s = cons_names[i];
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // check primitive variables
    std::map<std::string,int> prim_tags;
    for (int i=0; i<+PrimIdx::NUM; ++i) {
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
    Vector<Real> S(+ConsIdx::NUM), Q(+PrimIdx::NUM);

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

                for (int n=0; n<+ConsIdx::NUM; ++n) {
                    S[n] = src4(i,j,k,n);
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

    Array<Real,+ConsIdx::NUM> U;

    Real max_speed = std::numeric_limits<Real>::max();

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls
        const FArrayBox& vfrac = eb_data.volfrac[mfi];

        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();
#endif
        Array4<const Real> const& data4 = data.array(mfi);

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif

                    for (int n = 0; n<+ConsIdx::NUM; ++n) {
                        U[i] = data4(i,j,k,n);
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
