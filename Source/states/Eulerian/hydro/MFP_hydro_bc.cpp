#include "MFP_hydro.H"
#include "MFP_hydro_defs.H"

Vector<set_bc> HydroState::bc_set = {
    &set_scalar_bc,
    &set_x_vel_bc,
    &set_y_vel_bc,
    &set_z_vel_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc
};

std::map<std::string, int> HydroState::bc_names = {{"interior",  PhysBCType::interior},
                                                   {"inflow",    PhysBCType::inflow},
                                                   {"outflow",   PhysBCType::outflow},
                                                   {"symmetry",  PhysBCType::symmetry},
                                                   {"slipwall",  PhysBCType::slipwall},
                                                   {"noslipwall",PhysBCType::noslipwall}};

void HydroState::wall_eb(int eb_idx)
{

    std::pair<int, int>& idx = eb_bc_index[eb_idx];

    switch(idx.first) {
    case +HydroDef::WallIndex::HydroSlipWall:
        slip_wall_eb();
        break;
    case +HydroDef::WallIndex::HydroNoSlipWall:
        no_slip_wall_eb(idx.second);
        break;
    case +HydroDef::WallIndex::HydroDefined:
        defined_wall_eb(idx.second);
        break;
    default:
        Abort("How did we get here?");
    }
}

void HydroState::slip_wall_eb()
{
    return;
}

void HydroState::no_slip_wall_eb(int idx)
{
    NoSlipWallData& dat = no_slip_wall_eb_data[idx];

    return;
}

void HydroState::defined_wall_eb(int idx)
{
    Array<Real,+HydroDef::PrimIdx::NUM>& dat = defined_wall_eb_data[idx];

    return;
}
