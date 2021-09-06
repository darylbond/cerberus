#include "MFP_field.H"
#include "MFP_fillbc.H"

Vector<set_bc> FieldState::bc_set = {
    &set_x_D_bc,
    &set_y_D_bc,
    &set_z_D_bc,
    &set_x_B_bc,
    &set_y_B_bc,
    &set_z_B_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc,
};

Vector<std::string> FieldState::cons_names = {
    "x_D",
    "y_D",
    "z_D",
    "x_B",
    "y_B",
    "z_B",
    "phi",
    "psi",
    "mu",
    "ep",
};

Array<int,+FieldDef::VectorIdx::Cons> FieldState::vector_idx = {+FieldDef::FieldDef::ConsIdx::Bx, +FieldDef::FieldDef::ConsIdx::Dx};

std::map<std::string, int> FieldState::bc_names = {{"interior",  PhysBCType::interior},
                                                   {"inflow",    PhysBCType::inflow},
                                                   {"outflow",   PhysBCType::outflow},
                                                   {"symmetry",  PhysBCType::symmetry},
                                                   {"asymmetry",  4}};

std::string FieldState::tag = "field";
bool FieldState::registered = GetStateFactory().Register(FieldState::tag, StateBuilder<FieldState>);


FieldState::FieldState(){}

FieldState::FieldState(const sol::table &def)
{
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
}

FieldState::~FieldState(){}

void FieldState::variable_setup()
{

    ParmParse pgm("geometry");
    Vector<int> is_per(AMREX_SPACEDIM,0);
    pgm.queryarr("is_periodic",is_per);

    Vector<std::string> comp_names(+FieldDef::ConsIdx::NUM);

    // get the list of bc setting functions that applies to this type of state

    boundary_conditions.fill_bc.resize(+FieldDef::ConsIdx::NUM);

    for (int icomp=0; icomp < +FieldDef::ConsIdx::NUM; ++icomp) {
        set_bc s = bc_set[icomp]; // the function that sets the bc

        // make sure our periodicity isn't being overwritten
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            if (is_per[d]) {
                boundary_conditions.phys_fill_bc[icomp].setLo(d, PhysBCType::interior);
                boundary_conditions.phys_fill_bc[icomp].setHi(d, PhysBCType::interior);
            }
        }

        // grab the per component BCRec and apply the set bc function to it
        (*s)(boundary_conditions.fill_bc[icomp], boundary_conditions.phys_fill_bc[icomp]);
    }

    for (int icomp=0; icomp < +FieldDef::ConsIdx::NUM; ++icomp) {
        comp_names[icomp] = cons_names[icomp] + "-" + name;
    }

    int ng = get_num_grow();

#ifdef AMREX_USE_EB
    Interpolater* interp = &eb_cell_cons_interp;
#else
    Interpolater* interp = &cell_cons_interp;
#endif

    data_idx = MFP::get_desc_lst().size();

    MFP::get_desc_lst().addDescriptor(data_idx, IndexType::TheCellType(),
                                      StateDescriptor::Point, ng, +FieldDef::ConsIdx::NUM,
                                      interp, false, true);

    MFP::get_desc_lst().setComponent(data_idx, 0, comp_names, boundary_conditions.fill_bc, FillBC());


    if (MFP::verbosity >= 1) {
        Print() << str();
    }
}

void FieldState::init_data(MFP* mfp)
{

    const Real* dx = mfp->Geom().CellSize();
    const Real* prob_lo = mfp->Geom().ProbLo();

    MultiFab& S_new = mfp->get_new_data(data_idx);

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        FArrayBox& src = S_new[mfi];

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);
        Array4<Real> const& h4 = src.array();

        Real x, y, z;
        for     (int k = lo.z; k <= hi.z; ++k) {
            z = prob_lo[2] + (k + 0.5)*dx[2];
            for   (int j = lo.y; j <= hi.y; ++j) {
                y = prob_lo[1] + (j + 0.5)*dx[1];
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    x = prob_lo[0] + (i + 0.5)*dx[0];

                    // grab the primitive variables as defined by the user functions
                    for (int icomp=0; icomp<+FieldDef::ConsIdx::NUM; ++icomp) {
                        const std::string& cons_name = cons_names[icomp];
                        const auto& f = functions[cons_name];

                        h4(i,j,k,icomp) = f(x, y, z);

                    }
                }
            }
        }
    }

    return;
}

Real FieldState::get_allowed_time_step(MFP* mfp) const
{
    const Real* dx = mfp->Geom().CellSize();

    Real dt = dx[0]/fastest_speed;
    for (int i=1; i<AMREX_SPACEDIM; ++i) {
        dt = std::min(dt, dx[i]/fastest_speed);
    }

    return dt;
}

void FieldState::get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                                 #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                                 #endif
                                 ) const
{
    BL_PROFILE("FieldState::get_state_values");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
    Array4<const Real> const& vf4 = vfrac.array();
#endif

    updated.resize(0);

    // check conserved variables
    std::map<std::string,int> cons_tags;
    for (int i=0; i<+FieldDef::ConsIdx::NUM; ++i) {
        const std::string s = cons_names[i];
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // additional variables
    Vector<std::string> other;

#ifdef AMREX_USE_EB
    const std::string vfrac_name = "vfrac-"+name;
    bool load_vfrac = out.find(vfrac_name) != out.end();
    if (load_vfrac) other.push_back(vfrac_name);
#endif

    updated.insert(updated.end(), other.begin(), other.end());

    std::map<std::string,Array4<Real>> out4;
    for (const std::string& s : updated) {
        out[s].resize(box, 1);
        out4[s] = out[s].array();
    }

    Array4<const Real> const& src4 = src.array();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vf4(i,j,k) == 0.0) {
                    for (const std::string& s : updated) {
                        out4[s](i,j,k) = 0.0;
                    }
                    continue;
                }
#endif

                if (!cons_tags.empty()) {
                    for (const auto& var : cons_tags) {
                        out4[var.first](i,j,k) = src4(i,j,k,var.second);
                    }
                }

#ifdef AMREX_USE_EB
                if (load_vfrac)  out4[vfrac_name](i,j,k)  = vf4(i,j,k);
#endif
            }
        }
    }


    return;
}

void FieldState::set_udf()
{
    sol::state& lua = MFP::lua;

    bool success;

    sol::table state_def = lua["states"][name];

    // check if we have 'value' defined
    sol::table value = state_def["value"].get_or(sol::table());

    if ((!value.valid() || value.empty()) && ParallelDescriptor::IOProcessor())
        Warning("WARNING: State '"+name+"' does not have 'value' defined for initial conditions, using defaults");

    // get a list of any initialisation functions that need to be called during the run

    sol::table dynamic = state_def["dynamic"].get_or(sol::table());

    const Vector<std::pair<int,Real>> init_with_value = {
        {+FieldDef::ConsIdx::phi, 0.0},
        {+FieldDef::ConsIdx::psi, 0.0},
        {+FieldDef::ConsIdx::mu, 1.0},
        {+FieldDef::ConsIdx::ep, 1.0},
    };

    for (int i = 0; i<cons_names.size(); ++i) {

        const std::string &comp = cons_names[i];

        Optional3D1VFunction v;

        if (!value.valid() || value.empty()) {
            v.set_value(0.0);
            success = false;
        } else {
            success = get_udf(value[comp], v, 0.0);
        }

        // set default values
        if (!success) {
            for (const auto &j : init_with_value) {
                if (i == j.first) {
                    v.set_value(j.second);
                    break;
                }
            }
        }

        functions[comp] = v;

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

void FieldState::init_from_lua()
{
    BL_PROFILE("FieldState::init_from_lua");
    sol::state& lua = MFP::lua;

    const sol::table state_def = lua["states"][name];

    is_static = (bool) state_def["static"].get_or(0);

    set_udf();

    EulerianState::init_from_lua();

    //
    // boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};

    // mapping between name and index for various groupings
    std::map<std::string,std::map<std::string, int>> bc2index;
    for (int i=+FieldDef::ConsIdx::Dx; i<=+FieldDef::ConsIdx::Dz; ++i) {
        bc2index["fill_D_bc"][cons_names[i]] = i;
    }

    for (int i=+FieldDef::ConsIdx::Bx; i<=+FieldDef::ConsIdx::Bz; ++i) {
        bc2index["fill_B_bc"][cons_names[i]] = i;
    }

    bc2index["fill_ep_bc"][cons_names[+FieldDef::ConsIdx::ep]] = +FieldDef::ConsIdx::ep;
    bc2index["fill_mu_bc"][cons_names[+FieldDef::ConsIdx::mu]] = +FieldDef::ConsIdx::mu;

    bc2index["fill_psi_bc"] = {{"psi",+FieldDef::ConsIdx::psi}};
    bc2index["fill_phi_bc"] = {{"phi",+FieldDef::ConsIdx::phi}};

    BoundaryState &bs = boundary_conditions;


    bs.fill_bc.resize(+FieldDef::ConsIdx::NUM);
    bs.phys_fill_bc.resize(+FieldDef::ConsIdx::NUM);

    for (int ax = 0; ax < AMREX_SPACEDIM; ++ax) {


        for (int lh=0; lh<2; ++lh) {

#ifdef AMREX_USE_EB
            bool is_symmetry = false;
#endif
            for (const auto &bc : bc2index) {

                // get the base boundary condition for cell centered values
                std::string side_bc = state_def["bc"][dir_name[ax]][side_name[lh]][bc.first].get_or<std::string>("outflow");
                int i_side_bc = bc_names.at(side_bc);
#ifdef AMREX_USE_EB
                if (i_side_bc == PhysBCType::symmetry || i_side_bc == 4) is_symmetry = true;
#endif
                // fill in the bc list for AMReX as well as gather any custom values/functions
                for (const auto &var : bc.second) {

                    if (lh==0) {
                        bs.phys_fill_bc[var.second].setLo(ax,i_side_bc);
                    } else {
                        bs.phys_fill_bc[var.second].setHi(ax,i_side_bc);
                    }

                    const sol::object bcv = state_def["bc"][dir_name[ax]][side_name[lh]][var.first].get_or(sol::object());

                    Optional3D1VFunction v;

                    // special case for phi and psi (set to zero in boundary)
                    if (var.second == +FieldDef::ConsIdx::phi || var.second == +FieldDef::ConsIdx::psi) {
                        get_udf(bcv,v,0.0);
                    } else {
                        v = get_udf(bcv);
                    }
                    bs.set(ax,cons_names[var.second],lh,v);

                    // special case for inflow condition
                    if (i_side_bc == PhysBCType::inflow && !v.is_valid()) {
                        Abort("Setting '"+bc.first+" = inflow' requires all primitive variables to be defined, '" + var.first + "' is not defined");
                    }
                }
            }

#ifdef AMREX_USE_EB
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

    relative_div_speed = state_def["div_transport"].get_or(0.0);

    div_speed = relative_div_speed*MFP::lightspeed;

    fastest_speed = std::max(MFP::lightspeed, div_speed);

    // or we can use the projection method for divergence error control
    project_divergence = state_def["project_divergence"].get_or(0);


    return;
}
