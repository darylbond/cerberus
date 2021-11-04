#include "MFP_mhd.H"
#include "MFP_lua.H"
#include "MFP.H"
#include "MFP_diagnostics.H"
#include "MFP_transforms.H"
#include "MFP_mhd_refine.H"

#include "Eigen"

Vector<std::string> MHDState::cons_names = {
    "rho",
    "x_mom",
    "y_mom",
    "z_mom",
    "nrg",
    "x_B",
    "y_B",
    "z_B",
    "psi"
};

Vector<std::string> MHDState::prim_names = {
    "rho",
    "x_vel",
    "y_vel",
    "z_vel",
    "p",
    "x_B",
    "y_B",
    "z_B",
    "psi"
};

Array<int,2> MHDState::cons_vector_idx = {+MHDDef::ConsIdx::Xmom, +MHDDef::ConsIdx::Bx};
Array<int,2> MHDState::prim_vector_idx = {+MHDDef::PrimIdx::Xvel, +MHDDef::PrimIdx::Bx};

std::map<std::string, int> MHDState::bc_names = {{"interior",  PhysBCType::interior},
                                                   {"inflow",    PhysBCType::inflow},
                                                   {"outflow",   PhysBCType::outflow},
                                                   {"symmetry",  PhysBCType::symmetry},
                                                   {"slipwall",  PhysBCType::slipwall},
                                                   {"noslipwall",PhysBCType::noslipwall}};

Vector<set_bc> MHDState::bc_set = {
    &set_scalar_bc,
    &set_x_vel_bc,
    &set_y_vel_bc,
    &set_z_vel_bc,
    &set_scalar_bc,
    &set_x_B_bc,
    &set_y_B_bc,
    &set_z_B_bc,
    &set_scalar_bc,
};

std::string MHDState::tag = "mhd";
bool MHDState::registered = GetStateFactory().Register(MHDState::tag, StateBuilder<MHDState>);

MHDState::MHDState(){}

MHDState::MHDState(const sol::table& def)
{
    num_grow = 0;
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
}

MHDState::~MHDState(){}

#ifdef AMREX_USE_EB
void MHDState::set_eb_bc(const sol::table &bc_def)
{

    std::string bc_type = bc_def.get<std::string>("type");

    if (bc_type == MHDSlipWall::tag) {
        eb_bcs.push_back(std::unique_ptr<MHDBoundaryEB>(new MHDSlipWall(flux_solver.get())));
    } else if (bc_type == DirichletWallMHD::tag) {
        eb_bcs.push_back(std::unique_ptr<MHDBoundaryEB>(new DirichletWallMHD(flux_solver.get(), bc_def)));
    } else {
        Abort("Requested EB bc of type '" + bc_type + "' which is not compatible with state '" + name + "'");
    }
}
#endif


void MHDState::set_udf()
{
    using namespace std::placeholders;

    sol::state& lua = MFP::lua;

    sol::table state_def = lua["states"][name];

    // check if we have 'value' defined
    const sol::table value = state_def["value"].get_or(sol::table());

    if (!value.valid())
        Abort("State "+name+" does not have 'value' defined for initial conditions");

    functions.resize(+MHDDef::PrimIdx::NUM);

    // now for the primitives
    for (int i = 0; i<+MHDDef::PrimIdx::NUM; ++i) {

        std::string comp = prim_names[i];

        Optional3D1VFunction v;

        get_udf(value[comp], v, 0.0);

        functions[i] = v;

    }

    return;
}

void MHDState::set_flux()
{

    if (!is_transported()) return;

    ClassFactory<MHDRiemannSolver> rfact = GetMHDRiemannSolverFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;
    state_def["gamma"] = gamma;
    state_def["div_transport"] = div_transport;

    std::string flux = state_def["flux"].get_or<std::string>("null");

    if (flux == "null")
        Abort("Flux option required for state '"+name+"'. Options are "+vec2str(rfact.getKeys()));

    flux_solver = rfact.Build(flux, state_def);

    if (!flux_solver)
        Abort("Invalid flux solver option '"+flux+"'. Options are "+vec2str(rfact.getKeys()));


    return;

}

void MHDState::set_shock_detector()
{

    ClassFactory<MHDShockDetector> sdfact = GetMHDShockDetectorFactory();

    sol::table sd_def = MFP::lua["states"][name]["shock_detector"].get_or(sol::table());

    if (!sd_def.valid()) return;

    sd_def["global_idx"] = global_idx;

    std::string sd_name = sd_def["name"].get_or<std::string>("");

    shock_detector = sdfact.Build(sd_name, sd_def);

    if (!sd_name.empty() && !shock_detector)
        Abort("Invalid shock_detector option '"+sd_name+"'. Options are "+vec2str(sdfact.getKeys()));
}

void MHDState::set_refinement()
{

    ClassFactory<Refinement> rfact = GetMHDRefinementFactory();

    sol::table r_def = MFP::lua["states"][name]["refinement"].get_or(sol::table());

    if (!r_def.valid()) return;

    r_def["global_idx"] = global_idx;

    std::string r_name = r_def["name"].get_or<std::string>("");

    refine = rfact.Build(r_name, r_def);

    if (!r_name.empty() && !refine)
        Abort("Invalid refinement option '"+r_name+"'. Options are "+vec2str(rfact.getKeys()));
}

void MHDState::init_from_lua()
{
    BL_PROFILE("MHDState::init_from_lua");

    EulerianState::init_from_lua();

    sol::state& lua = MFP::lua;

    const sol::table state_def = lua["states"][name];


    //
    // get gamma
    //


    gamma = state_def["gamma"];


    //
    // user defined functions
    //
    set_udf();


    //
    // domain boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};
    const Vector<std::string>& mhd_var = prim_names;
    const int N = mhd_var.size();

    BoundaryState &bs = boundary_conditions;
    bs.phys_fill_bc.resize(+MHDDef::PrimIdx::NUM);

    for (int ax = 0; ax < AMREX_SPACEDIM; ++ax) {
        for (int lh=0; lh<2; ++lh) {

            std::string side_bc = state_def["bc"][dir_name[ax]][side_name[lh]]["fill_mhd_bc"].get_or<std::string>("outflow");
            int i_side_bc = bc_names.at(side_bc);

            // get any custom values/functions
            for (int j=0; j<N; ++j) {

                if (lh==0) {
                    bs.phys_fill_bc[j].setLo(ax,i_side_bc);
                } else {
                    bs.phys_fill_bc[j].setHi(ax,i_side_bc);
                }

                const sol::object v = state_def["bc"][dir_name[ax]][side_name[lh]][mhd_var[j]].get_or(sol::object());
                Optional3D1VFunction f = get_udf(v);
                bs.set(ax,mhd_var[j],lh,f);

                // special case for inflow condition
                if (i_side_bc == PhysBCType::inflow && !f.is_valid()) {
                    Abort("Setting 'fill_mhd_bc = inflow' requires all primitive variables to be defined, '" + mhd_var[j] + "' is not defined");
                }
            }



#ifdef AMREX_USE_EB
            bool is_symmetry = (i_side_bc == PhysBCType::symmetry) || (i_side_bc == PhysBCType::slipwall);
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
    // divergence handling
    //

    div_transport = state_def["div_transport"].get_or(0.0);

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

void MHDState::variable_setup(Vector<int> periodic)
{

    boundary_conditions.fill_bc.resize(+MHDDef::PrimIdx::NUM);

    for (int icomp=0; icomp < +MHDDef::ConsIdx::NUM; ++icomp) {
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

    Vector<std::string> comp_names(+MHDDef::ConsIdx::NUM);
    for (int icomp=0; icomp < +MHDDef::ConsIdx::NUM; ++icomp) {
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
                           StateDescriptor::Point, ng, +MHDDef::ConsIdx::NUM,
                           interp, state_data_extrap,
                           store_in_checkpoint);

    desc_lst.setComponent(
                data_idx, 0, comp_names, boundary_conditions.fill_bc,
                FillBC());


    if (MFP::verbosity >= 1) {
        Print() << str();
    }
}

void MHDState::init_data(MFP* mfp, const Real time)
{

    const Real* dx = mfp->Geom().CellSize();
    const Real* prob_lo = mfp->Geom().ProbLo();

    MultiFab& S_new = mfp->get_data(data_idx, time);

    Array<Real,+MHDDef::ConsIdx::NUM> U;
    Array<Real,+MHDDef::PrimIdx::NUM> Q;

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        FArrayBox& S_data = S_new[mfi];

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);
        Array4<Real> const& S4 = S_data.array();

#ifdef AMREX_USE_EB
        FabArray<EBCellFlagFab>& flags = mfp->get_eb_data(global_idx).flags;
        Array4<const EBCellFlag> const& flag4 = flags.array(mfi);
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
                        for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
                            S4(i,j,k,n) = 0.0;
                        }
                        continue;
                    }
#endif

                    // grab the primitive variables as defined by the user functions
                    for (int icomp=0; icomp<+MHDDef::PrimIdx::NUM; ++icomp) {
                        const auto& f = functions[icomp];

                        Q[icomp] = f(x, y, z);

                    }

                    // convert primitive to conserved
                    prim2cons(Q, U);

                    // copy into array
                    for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
                        S4(i,j,k,n) = U[n];
                    }
                }
            }
        }
    }
}

// in place conversion from conserved to primitive
bool MHDState::cons2prim(Array<Real, +MHDDef::ConsIdx::NUM>& U, Array<Real, +MHDDef::PrimIdx::NUM>& Q) const
{
    BL_PROFILE("MHDState::cons2prim");

    Real rho = U[+MHDDef::ConsIdx::Density];
    Real rhoinv = 1/rho;
    Real u = U[+MHDDef::ConsIdx::Xmom]*rhoinv;
    Real v = U[+MHDDef::ConsIdx::Ymom]*rhoinv;
    Real w = U[+MHDDef::ConsIdx::Zmom]*rhoinv;

    Real Bx = U[+MHDDef::ConsIdx::Bx];
    Real By = U[+MHDDef::ConsIdx::By];
    Real Bz = U[+MHDDef::ConsIdx::Bz];

    Real nrg = U[+MHDDef::ConsIdx::Eden] - 0.5*(rho*(u*u + v*v + w*w) + Bx*Bx + By*By + Bz*Bz);

    Real p = nrg*(gamma - 1.0);

    Q[+MHDDef::PrimIdx::Density] = rho;
    Q[+MHDDef::PrimIdx::Xvel] = u;
    Q[+MHDDef::PrimIdx::Yvel] = v;
    Q[+MHDDef::PrimIdx::Zvel] = w;
    Q[+MHDDef::PrimIdx::Prs] = p;
    Q[+MHDDef::PrimIdx::Bx] = Bx;
    Q[+MHDDef::PrimIdx::By] = By;
    Q[+MHDDef::PrimIdx::Bz] = Bz;
    Q[+MHDDef::PrimIdx::psi] = U[+MHDDef::ConsIdx::psi];


    return prim_valid(Q);
}

void MHDState::prim2cons(Array<Real, +MHDDef::PrimIdx::NUM>& Q, Array<Real, +MHDDef::ConsIdx::NUM>& U) const
{
    BL_PROFILE("MHDState::prim2cons");

    Real rho = Q[+MHDDef::PrimIdx::Density];
    Real mx = rho*Q[+MHDDef::PrimIdx::Xvel];
    Real my = rho*Q[+MHDDef::PrimIdx::Yvel];
    Real mz = rho*Q[+MHDDef::PrimIdx::Zvel];
    Real p = Q[+MHDDef::PrimIdx::Prs];

    Real Bx = Q[+MHDDef::PrimIdx::Bx];
    Real By = Q[+MHDDef::PrimIdx::By];
    Real Bz = Q[+MHDDef::PrimIdx::Bz];

    U[+MHDDef::ConsIdx::Density] = rho;
    U[+MHDDef::ConsIdx::Xmom] = mx;
    U[+MHDDef::ConsIdx::Ymom] = my;
    U[+MHDDef::ConsIdx::Zmom] = mz;
    U[+MHDDef::ConsIdx::Eden] = p/(gamma - 1) + 0.5*(mx*mx + my*my + mz*mz)/rho + 0.5*(Bx*Bx + By*By + Bz*Bz);
    U[+MHDDef::ConsIdx::Bx] = Bx;
    U[+MHDDef::ConsIdx::By] = By;
    U[+MHDDef::ConsIdx::Bz] = Bz;
    U[+MHDDef::ConsIdx::psi] = Q[+MHDDef::PrimIdx::psi];

}

bool MHDState::prim_valid(const Array<Real, +MHDDef::PrimIdx::NUM> &Q) const
{
    if ((Q[+MHDDef::PrimIdx::Density] <= 0.0) ||  (Q[+MHDDef::PrimIdx::Prs] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

bool MHDState::cons_valid(const Array<Real, +MHDDef::ConsIdx::NUM> &U) const
{
    if ((U[+MHDDef::ConsIdx::Density] <= 0.0) ||  (U[+MHDDef::ConsIdx::Eden] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

RealArray MHDState::get_speed_from_cons(const Array<Real,+MHDDef::ConsIdx::NUM> &U) const
{
    BL_PROFILE("MHDState::get_speed_from_cons");
    Real rho = U[+MHDDef::ConsIdx::Density];
    Real ed = U[+MHDDef::ConsIdx::Eden];
    Real rhoinv = 1/rho;

    Real ux = U[+MHDDef::ConsIdx::Xmom]*rhoinv;
    Real uy = U[+MHDDef::ConsIdx::Ymom]*rhoinv;
    Real uz = U[+MHDDef::ConsIdx::Zmom]*rhoinv;

    Real kineng = 0.5*rho*(ux*ux + uy*uy + uz*uz);
    Real p = (ed - kineng)*(gamma - 1);

    Real Bx = U[+MHDDef::ConsIdx::Bx];
    Real By = U[+MHDDef::ConsIdx::By];
    Real Bz = U[+MHDDef::ConsIdx::Bz];

    Real B = std::sqrt(Bx*Bx + By*By + Bz*Bz);

    Real cf = std::sqrt( (gamma*p + B) / rho);

    RealArray s = {AMREX_D_DECL(cf + std::abs(ux), cf + std::abs(uy), cf + std::abs(uz))};

    return s;

}

RealArray MHDState::get_speed_from_prim(const Array<Real, +MHDDef::PrimIdx::NUM> &Q) const
{
    BL_PROFILE("MHDState::get_speed_from_prim");

    Real rho = Q[+MHDDef::PrimIdx::Density];
    Real ux = Q[+MHDDef::PrimIdx::Xvel];
    Real uy = Q[+MHDDef::PrimIdx::Yvel];
    Real uz = Q[+MHDDef::PrimIdx::Zvel];
    Real p = Q[+MHDDef::PrimIdx::Prs];

    Real Bx = Q[+MHDDef::PrimIdx::Bx];
    Real By = Q[+MHDDef::PrimIdx::By];
    Real Bz = Q[+MHDDef::PrimIdx::Bz];

    Real B = std::sqrt(Bx*Bx + By*By + Bz*Bz);

    Real cf = std::sqrt( (gamma*p + B) / rho);

    RealArray s = {AMREX_D_DECL(cf + std::abs(ux), cf + std::abs(uy), cf + std::abs(uz))};

    return s;

}

Vector<std::string> MHDState::get_plot_output_names() const
{
    Vector<std::string> out;
    out.insert(out.end(), cons_names.begin(), cons_names.end());
    out.insert(out.end(), prim_names.begin(), prim_names.end());
#ifdef AMREX_USE_EB
    out.push_back("vfrac");
#endif

    return out;
}

void MHDState::get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                                 #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                                 #endif
                                 ) const
{
    BL_PROFILE("MHDState::get_state_values");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
    Array4<const Real> const& vf4 = vfrac.array();
#endif

    updated.resize(0);

    // check conserved variables
    std::map<std::string,int> cons_tags;
    for (int i=0; i<+MHDDef::ConsIdx::NUM; ++i) {
        const std::string s = cons_names[i];
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // check primitive variables
    std::map<std::string,int> prim_tags;
    for (int i=0; i<+MHDDef::PrimIdx::NUM; ++i) {
        const std::string s = prim_names[i];
        if (s == cons_names[0]) continue;
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        prim_tags[var_name] = i;
        updated.push_back(var_name);
    }

#ifdef AMREX_USE_EB
    const std::string vfrac_name = "vfrac-"+name;
    bool load_vfrac = out.find(vfrac_name) != out.end();
#endif

    std::map<std::string,Array4<Real>> out4;
    for (const std::string& s : updated) {
        out[s].resize(box, 1);
        out[s].setVal(0.0);
        out4[s] = out[s].array();
    }

    // temporary storage for retrieving the state data
    Array<Real,+MHDDef::ConsIdx::NUM> S;
    Array<Real,+MHDDef::PrimIdx::NUM> Q;

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

                for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
                    S[n] = src4(i,j,k,n);

                    if (std::isnan(S[n])) {
                        Abort();
                    }

                }


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

Real MHDState::get_allowed_time_step(MFP* mfp) const
{
    const Real* dx = mfp->Geom().CellSize();

    MultiFab& data = mfp->get_new_data(data_idx);

    Array<Real,+MHDDef::ConsIdx::NUM> U;

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

                    for (int n = 0; n<+MHDDef::ConsIdx::NUM; ++n) {
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

    if (div_transport > 0.0) {
        ParallelDescriptor::ReduceRealMin(dt);
        flux_solver->update_max_speed(mfp->cfl*dx[0]/dt);
    }

    return dt;
}

void MHDState::calc_velocity(const Box& box,
                               FArrayBox& cons,
                               FArrayBox &prim
                               #ifdef AMREX_USE_EB
                               ,const FArrayBox& vfrac
                               #endif
                               ) const
{
    BL_PROFILE("MHDState::calc_velocity");


    prim.resize(box, AMREX_SPACEDIM);

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);
    Array4<Real> const& s4 = cons.array();
    Array4<Real> const& p4 = prim.array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();
#endif

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vfrac4(i,j,k) == 0.0) continue;
#endif

                const Real irho = 1/s4(i,j,k,+MHDDef::ConsIdx::Density);

                for (int n = 0; n<AMREX_SPACEDIM; ++n) {
                    p4(i,j,k,n) =  irho*s4(i,j,k,+MHDDef::ConsIdx::Xmom+n);
                }

            }
        }
    }

    return;
}

void MHDState::calc_primitives(const Box& box,
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
    BL_PROFILE("MHDState::calc_primitives");

    Array<Real,+MHDDef::ConsIdx::NUM> U;
    Array<Real,+MHDDef::PrimIdx::NUM> Q;

    prim.resize(box, +MHDDef::PrimIdx::NUM);

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
                            for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
                                U[n] += vf*s4(ii,jj,kk,n);
                            }

                            vtot += vf;
                        }
                    }

                    // if we were close enough to a boundary to have valid data we
                    // average out the volume fraction weighted contributions, otherwise,
                    // fill in the primitives with zeros
                    if (vtot > 0.0) {
                        for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
                            U[n] /= vtot;
                        }
                    } else {
                        for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
                            p4(i,j,k,n) = 0.0;
                        }
                        continue;
                    }
                } else {
#endif
                    // grab the conserved variables
                    for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
                        U[n] = s4(i,j,k,n);
                    }
#ifdef AMREX_USE_EB
                }
#endif


                // convert to primitive
                cons2prim(U, Q);

                // copy into primitive
                for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
                    p4(i,j,k,n) = Q[n];
                }
            }
        }
    }

    return;
}

void MHDState::update_boundary_cells(const Box& box,
                                       const Geometry &geom,
                                       FArrayBox &prim,
                                       #ifdef AMREX_USE_EB
                                       const FArrayBox& vfrac,
                                       #endif
                                       const Real time) const
{
    BL_PROFILE("MHDState::update_boundary_cells");
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
                        for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
                            Q[prim_names[n]] = p4(grab[0],grab[1],grab[2],n);
                        }

                        // update the primitives from our UDFs, but only those that are valid functions
                        for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
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

void MHDState::calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                     #ifdef AMREX_USE_EB
                                     ,const EBCellFlagFab &flag
                                     ,const FArrayBox &vfrac
                                     #endif
                                     ) const
{
    BL_PROFILE("MHDState::calc_reconstruction");

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& src4 = prim.array();

#ifdef AMREX_USE_EB

    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);


    Array4<const EBCellFlag> const& f4 = flag.array();
    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    Vector<Real> stencil(reconstructor->stencil_length);
    int offset = reconstructor->stencil_length/2;
    Array<int,3> stencil_index;
    Real lo_face, hi_face;


    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        // make sure our arrays for putting lo and hi reconstructed values into
        // are the corect size
        rlo[d].resize(box, +MHDDef::ConsIdx::NUM);
        rhi[d].resize(box, +MHDDef::ConsIdx::NUM);

#ifdef AMREX_USE_EB
        if (check_eb) {
            rlo[d].copy(prim,box);
            rhi[d].copy(prim,box);
        }
#endif

        Array4<Real> const& lo4 = rlo[d].array();
        Array4<Real> const& hi4 = rhi[d].array();

        // cycle over all components
        for (int n = 0; n<+MHDDef::ConsIdx::NUM; ++n) {

            for     (int k = lo.z; k <= hi.z; ++k) {
                for   (int j = lo.y; j <= hi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                            for (int i = lo.x; i <= hi.x; ++i) {



#ifdef AMREX_USE_EB
                        if (check_eb) {
                            // covered cell doesn't need calculating
                            if (f4(i,j,k).isCovered() || check_covered_stencil(f4, i, j, k, d, reconstructor->stencil_length)) {
                                continue;
                            }
                        }
#endif
                        stencil_index.fill(0);
                        for (int s=0; s<reconstructor->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], n);
                        }

                        // perform reconstruction
                        reconstructor->get_face_values(stencil, lo_face, hi_face);

                        lo4(i,j,k,n) = lo_face;
                        hi4(i,j,k,n) = hi_face;
                    }
                }
            }
        }
    }


    return;
}

void MHDState::calc_time_averaged_faces(const Box& box,
                                          const FArrayBox &prim,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                          #ifdef AMREX_USE_EB
                                          const EBCellFlagFab& flag,
                                          #endif
                                          const Real* dx,
                                          Real dt) const
{
    BL_PROFILE("MHDState::calc_time_averaged_faces");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& p4 = prim.array();

    Array<Real,+MHDDef::PrimIdx::NUM> Q;

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

                for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
                    Q[n] = p4(i,j,k,n);
                }

                const RealArray local_c = get_speed_from_prim(Q);


                for (int d=0; d<AMREX_SPACEDIM; ++d) {

                    Array4<Real> const& lo4 = rlo[d].array();
                    Array4<Real> const& hi4 = rhi[d].array();

                    dt_2dx = 0.5*dt/dx[d];

                    // cycle over all components
                    for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n) {

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

void MHDState::face_bc(const int dir,
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
    BL_PROFILE("MHDState::face_bc");
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
        for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
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
                        for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
                            Q[prim_names[n]] = s4(i,j,k,n);
                        }

                        // update the primitives from our UDFs, but only those that are valid functions
                        for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
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
void MHDState::load_state_for_flux(const Array4<const Real> &face,
                                     int i, int j, int k,
                                   Array<Real,+MHDDef::PrimIdx::NUM> &S) const
{
    BL_PROFILE("MHDState::load_state_for_flux");

    for (size_t n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
        S[n] = face(i,j,k,n);
    }

    return;
}

void MHDState::calc_fluxes(const Box& box,
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
    BL_PROFILE("MHDState::calc_fluxes");

    Array<int, 3> index;
    Array<Real,+MHDDef::PrimIdx::NUM> L, R;
    Array<Real,+MHDDef::ConsIdx::NUM> F;

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

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
                    transform_global2local(L, d, prim_vector_idx);
                    transform_global2local(R, d, prim_vector_idx);


                    Real shk = 0.0;

                    if (shock_detector) {
                        shk = shock_detector->solve(L, R);
                    }

                    flux_solver->solve(L, R, F, &shk);

                    // rotate the flux back to local frame
                    transform_local2global(F, d, cons_vector_idx);

                    // load the flux into the array
                    for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
                        flux4(i,j,k,n) += F[n];

                        AMREX_ASSERT(std::isfinite(flux4(i,j,k,n)));

                    }

                }
            }
        }
    }

    return;
}


void MHDState::correct_face_prim(const Box& box,
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

    Array<Real,+MHDDef::PrimIdx::NUM> L_prim, R_prim;
    Array<Real,+MHDDef::ConsIdx::NUM> L_cons, R_cons;

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
                        for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n ) {
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
                        for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n ) {
                            loF = flux4(i-idx1[0],j-idx1[1],k-idx1[2],n);
                            hiF = flux4(i-idx1[0]+idx2[0],j-idx1[1]+idx2[1],k-idx1[2]+idx2[2],n);

                            L_cons[n] += 0.5*dt/dx[d2]*(loF - hiF);
                        }

                        // right side
                        for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n ) {
                            loF = flux4(i,j,k,n);
                            hiF = flux4(i+idx2[0],j+idx2[1],k+idx2[2],n);

                            R_cons[n] += 0.5*dt/dx[d2]*(loF - hiF);
                        }

                        // convert to primitive
                        cons2prim(L_cons, L_prim);
                        cons2prim(R_cons, R_prim);

                        // update reconstruction
                        for (int n=0; n<+MHDDef::PrimIdx::NUM; ++n ) {
                            hi4(i-idx1[0],j-idx1[1],k-idx1[2],n) = L_prim[n];
                            lo4(i,j,k,n) = R_prim[n];
                        }
                    }
                }
            }
        }
    }
}


#ifdef AMREX_USE_EB

void MHDState::calc_wall_fluxes(const Box& box,
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
    BL_PROFILE("MHDState::calc_wall_fluxes");

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& p4 = prim.array();
    Array4<const Real> const& bcent4 = bcent.array();
    Array4<const Real> const& bnorm4 = bnorm.array();

    Array<Array4<Real>,AMREX_SPACEDIM> flux4;
    Array<Array4<const Real>,AMREX_SPACEDIM> afrac4;

    Array<Vector<Real>,AMREX_SPACEDIM> wall_flux;
    for (int i=0; i<AMREX_SPACEDIM;++i) { wall_flux[i].resize(+MHDDef::ConsIdx::NUM);}

    Array4<const EBCellFlag> const& flag4 = flag.array();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
        afrac4[d] = afrac[d]->array();

        // zero out the flux accumulator
        fluxes[d].setVal(0.0);
    }

    const Array4<const Real>& bc_idx4 = bc_idx.array();

    Array<Real,+MHDDef::PrimIdx::NUM> cell_state;

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
                    const MHDBoundaryEB& bc = *eb_bcs[ebi];

                    // calculate the wall flux
                    bc.solve(wall_coord, wall_centre, cell_state, p4, i, j, k, dx, wall_flux);

                    // load the flux into the fab
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
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


void MHDState::write_info(nlohmann::json &js) const
{

    EulerianState::write_info(js);

    // write out stuff that is common to all states

    js["type"] = tag;
    js["type_idx"] = +get_type();

    js["gamma"] = gamma;

}