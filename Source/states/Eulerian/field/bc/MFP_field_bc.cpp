#ifdef AMREX_USE_EB

#include "MFP_field_bc.H"
#include "MFP_field.H"
#include "MFP_transforms.H"
#include "MFP_utility.H"

//-----------------------------------------------------------------------------

void FieldState::set_eb_bc(const sol::table &bc_def)
{
    std::string bc_type = bc_def.get_or<std::string>("type", CollectionWall::tag);

    if (bc_type == ConductingWall::tag) {
        eb_bcs.push_back(std::unique_ptr<FieldBoundaryEB>(new ConductingWall(flux_solver.get(), bc_def)));
    } else if (bc_type == ScalarPotentialWall::tag) {
        eb_bcs.push_back(std::unique_ptr<FieldBoundaryEB>(new ScalarPotentialWall(bc_def)));
    } else if (bc_type == SurfaceChargeWall::tag) {
        eb_bcs.push_back(std::unique_ptr<FieldBoundaryEB>(new SurfaceChargeWall(bc_def)));
    } else if (bc_type == CollectionWall::tag) {
        eb_bcs.push_back(std::unique_ptr<FieldBoundaryEB>(new CollectionWall(bc_def)));
    } else if (bc_type == DefinedWall::tag) {
        eb_bcs.push_back(std::unique_ptr<FieldBoundaryEB>(new DefinedWall(flux_solver.get(), bc_def)));
    } else {
        Abort("Requested EB bc of type '" + bc_type + "' which is not compatible with state '" + name + "'");
    }
}

//-----------------------------------------------------------------------------

FieldBoundaryEB::FieldBoundaryEB(){}
FieldBoundaryEB::~FieldBoundaryEB(){}

//-----------------------------------------------------------------------------

std::string DefinedWall::tag = "defined";

DefinedWall::DefinedWall(){}
DefinedWall::~DefinedWall(){}

DefinedWall::DefinedWall(FieldRiemannSolver *flux, const sol::table &bc_def)
{

    flux_solver = flux;

    // grab the wall state from the lua definition
    // only get stuff that is defined, anything else will
    // be filled in from the fluid side of the wall
    int i = 0;
    for (const auto &name : FieldState::cons_names) {
        sol::object val = bc_def[name];
        if (val.valid()) {
            wall_value.push_back({i, val.as<Real>()});
        }
        ++i;
    }


}

void DefinedWall::solve(Array<Array<Real,3>,3> &wall_coord,
                        Array<Real,+FieldDef::ConsIdx::NUM> &state,
                        Array<Array<Real,+FieldDef::ConsIdx::NUM>,AMREX_SPACEDIM> &F,
                        const Real* dx) const
{

    transform_global2local(state, wall_coord,  FieldState::vector_idx);

    // fabricate a state for inside the wall based on the provided state
    Array<Real, +FieldDef::ConsIdx::NUM> W = state;

    for (const auto& pair : wall_value) {
        W[pair.first] = pair.second;
    }

    Array<Real,+FieldDef::ConsIdx::NUM> normal_flux;
    flux_solver->solve(state, W, normal_flux);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, FieldState::vector_idx);

    // split it up into components
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for (int n=0; n<normal_flux.size(); ++n) {
            F[d][n] = wall_coord[0][d]*normal_flux[n];
        }
    }

}




//-----------------------------------------------------------------------------

std::string ConductingWall::tag = "conductor";

ConductingWall::ConductingWall(){}
ConductingWall::~ConductingWall(){}

ConductingWall::ConductingWall(FieldRiemannSolver* flux, const sol::table &bc_def)
{

    flux_solver = flux;

    // grab any specified normal and tangential fields
    B1_defined = false;
    sol::optional<Real> B_t1 = bc_def["B_t1"];
    if (B_t1) {
        B1_defined = true;
        wall_B1 = B_t1.value();
    } else {
        wall_B1 = 0.0;
    }

    B2_defined = false;
    sol::optional<Real> B_t2 = bc_def["B_t2"];
    if (B_t2) {
        B2_defined = true;
        wall_B2 = B_t2.value();
    } else {
        wall_B2 = 0.0;
    }

    D_defined = false;
    sol::optional<Real> D_n = bc_def["D_n"];
    if (D_n) {
        D_defined = true;
        wall_D = D_n.value();
    } else {
        wall_D = 0.0;
    }
}

void ConductingWall::solve(Array<Array<Real,3>,3> &wall_coord,
                           Array<Real,+FieldDef::ConsIdx::NUM> &state,
                           Array<Array<Real,+FieldDef::ConsIdx::NUM>,AMREX_SPACEDIM> &F,
                           const Real* dx) const
{

    //
    // get the inviscid flux
    //

    transform_global2local(state, wall_coord,  FieldState::vector_idx);

    // fabricate a state for inside the wall based on the provided state
    Array<Real, +FieldDef::ConsIdx::NUM> W = state;

    // https://en.wikipedia.org/wiki/Interface_conditions_for_electromagnetic_fields

    //(B2 - B1).n12 = 0, where B2 & B1 are the B vectors and n12 is the normal from 1->2
    // but B2=0 so B1.n12=0 thus BxR=-BxL
    W[+FieldDef::ConsIdx::Bx] *= -1;
    W[+FieldDef::ConsIdx::psi] = 0;

    if (B1_defined)
        W[+FieldDef::ConsIdx::By] = wall_B1;
    if (B2_defined)
        W[+FieldDef::ConsIdx::Bz] = wall_B2;

    // need to implement surface charge effects
    //(D2 - D1).n12 = sc, where D2 & D1 are the D vectors and n12 is the normal from 1->2
    // but D2=0 so D1.n12=-sc thus DxR=-DxL+sc
    // we assume here that we have a surface charge sufficient to allow for DxR=DxL
    W[+FieldDef::ConsIdx::Dy] *= -1;
    W[+FieldDef::ConsIdx::Dz] *= -1;
    W[+FieldDef::ConsIdx::phi] = 0;

    if (D_defined)
        W[+FieldDef::ConsIdx::Dx] = wall_D;

    Array<Real, +FieldDef::ConsIdx::NUM> normal_flux;
    flux_solver->solve(state, W, normal_flux);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, FieldState::vector_idx);

    // split it up into components
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for (int n=0; n<normal_flux.size(); ++n) {
            F[d][n] = wall_coord[0][d]*normal_flux[n];
        }
    }

}



//-----------------------------------------------------------------------------

std::string ScalarPotentialWall::tag = "scalar_potential";

ScalarPotentialWall::ScalarPotentialWall(){}
ScalarPotentialWall::~ScalarPotentialWall(){}

ScalarPotentialWall::ScalarPotentialWall(const sol::table &bc_def)
{
    get_udf(bc_def["phi"], phi, 0.0);
}

Vector<Vector<Real> > ScalarPotentialWall::get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const
{

    // calculated at face centre

    Vector<Vector<Real>> out(+EBType::NUM);
    out[+EBType::ScalarPotential] = {phi(0,0,0,t)};

    return out;
}

//-----------------------------------------------------------------------------

std::string SurfaceChargeWall::tag = "surface_charge";

SurfaceChargeWall::SurfaceChargeWall(){}
SurfaceChargeWall::~SurfaceChargeWall(){}

SurfaceChargeWall::SurfaceChargeWall(const sol::table &bc_def)
{
    get_udf(bc_def["charge"], charge_density, 0.0);
}

Vector<Vector<Real>> SurfaceChargeWall::get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const
{
    // needs to be calculated for cell centre

    Vector<Vector<Real>> out(+EBType::NUM);
    out[+EBType::SurfaceCharge] = {charge_density(0,0,0,t)};

    return out;
}

//-----------------------------------------------------------------------------

std::string SurfaceCurrentWall::tag = "surface_current";

SurfaceCurrentWall::SurfaceCurrentWall(){}
SurfaceCurrentWall::~SurfaceCurrentWall(){}

SurfaceCurrentWall::SurfaceCurrentWall(const sol::table &bc_def)
{
    get_udf(bc_def["j1"], current_1, 0.0);
    get_udf(bc_def["j2"], current_2, 0.0);
}

Vector<Vector<Real>> SurfaceCurrentWall::get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const
{

    Array<Real,3> J = {0.0, current_1(0,0,0,t), current_2(0,0,0,t)};

    Array<int,1> index = {0};

    transform_local2global(J, wall_coord, index);

    // needs to be calculated for cell centre

    Vector<Vector<Real>> out(+EBType::NUM);
    out[+EBType::SurfaceCurrent] = arr2vec(J);

    return out;
}

//-----------------------------------------------------------------------------

std::string VectorPotentialWall::tag = "vector_potential";

VectorPotentialWall::VectorPotentialWall(){}
VectorPotentialWall::~VectorPotentialWall(){}

VectorPotentialWall::VectorPotentialWall(const sol::table &bc_def)
{
    align_with_boundary = bc_def.get_or("align_with_boundary", false);
    get_udf(bc_def["A0"], A_0, 0.0);
    get_udf(bc_def["A1"], A_1, 0.0);
    get_udf(bc_def["A2"], A_2, 0.0);

    has_function = A_0.has_func() || A_1.has_func() || A_2.has_func();
}

Vector<Vector<Real>> VectorPotentialWall::get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const
{

    std::map<std::string, Real> data;

    // only load data if the boundary conditions are functions
    if (has_function) {

        AMREX_D_TERM(data["x"] = wall_centre[0];,data["y"] = wall_centre[1];,data["z"] = wall_centre[2];)

        data["nx"] = wall_coord[0][0];
        data["ny"] = wall_coord[0][1];
        data["nz"] = wall_coord[0][2];

        data["t1x"] = wall_coord[1][0];
        data["t1y"] = wall_coord[1][1];
        data["t1z"] = wall_coord[1][2];

        data["t2x"] = wall_coord[1][0];
        data["t2y"] = wall_coord[1][1];
        data["t2z"] = wall_coord[1][2];

        data["t"] = t;
    }

    Array<Real,3> A = {A_0(data), A_1(data), A_2(data)};

    if (align_with_boundary) {
        Array<int,1> index = {0};
        transform_local2global(A, wall_coord, index);
    }

    // calculated at face centre
    Vector<Vector<Real>> out(+EBType::NUM);
    out[+EBType::VectorPotential] = arr2vec(A);

    return out;
}

//-----------------------------------------------------------------------------

std::string CollectionWall::tag = "collection";

CollectionWall::CollectionWall(){}
CollectionWall::~CollectionWall(){}

CollectionWall::CollectionWall(const sol::table &bc_def)
{
    sol::table types = bc_def["types"];

    for (const auto& type : types) {
        std::string bc_type = type.second.as<std::string>();

        if (bc_type == ScalarPotentialWall::tag) {
            bcs.push_back(std::unique_ptr<FieldBoundaryEB>(new ScalarPotentialWall(bc_def)));
        } else if (bc_type == VectorPotentialWall::tag) {
            bcs.push_back(std::unique_ptr<FieldBoundaryEB>(new VectorPotentialWall(bc_def)));
        } else if (bc_type == SurfaceChargeWall::tag) {
            bcs.push_back(std::unique_ptr<FieldBoundaryEB>(new SurfaceChargeWall(bc_def)));
        } else if (bc_type == SurfaceCurrentWall::tag) {
            bcs.push_back(std::unique_ptr<FieldBoundaryEB>(new SurfaceCurrentWall(bc_def)));
        } else {
            Abort("Requested EB bc of type '" + bc_type + "' which is not compatible with 'collection'");
        }
    }
}

Vector<Vector<Real>> CollectionWall::get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const
{

    Vector<Vector<Real>> out(+EBType::NUM);

    for (const auto &bc : bcs) {
        const auto val = bc->get_wall_state(wall_centre, wall_coord, t);
        out[+bc->get_type()] = val[+bc->get_type()];
    }

    return out;
}

const bool CollectionWall::time_varying() const
{
    for (const auto &bc : bcs) {
        if (bc->time_varying()) {
            return true;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------


#endif
