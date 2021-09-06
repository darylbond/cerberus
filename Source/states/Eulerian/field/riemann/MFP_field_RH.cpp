#include "MFP_field_RH.H"
#include "MFP_utility.H"
#include "MFP.H"
#include "MFP_field.H"

//================================================================================

// J. Moreno, E. Oliva, P. Velarde, J.C.P. 2020, In Press

std::string FieldRH::tag = "RankineHugoniot";
bool FieldRH::registered = GetRiemannSolverFactory().Register(FieldRH::tag, RiemannSolverBuilder<FieldRH>);

FieldRH::FieldRH(){}
FieldRH::FieldRH(const int i)
{
    idx = i;

    FieldState& istate = FieldState::get_state(idx);

    c0 = MFP::lightspeed;
    ch = istate.div_speed;
    ch2 = ch*ch;
}

void FieldRH::solve(Array<Real,+FieldState::ConsIdx::NUM> &L,
                    Array<Real,+FieldState::ConsIdx::NUM> &R,
                    Array<Real,+FieldState::ConsIdx::NUM> &F) const
{
    BL_PROFILE("FieldRH::solve");
    std::fill(F.begin(), F.end(), 0);

    Real Dr, Dl; // y- & z- D fields
    Real Br, Bl; // y- & z- B fields
    Real Pl, Pr; // div correction factors
    Real fl, fr; // x- fields
    Real epr, epl;
    Real mur, mul;
    Real cr, cl;

    // D-wave

    Dl = L[+FieldState::ConsIdx::Dy];
    Dr = R[+FieldState::ConsIdx::Dy];

    epr = R[+FieldState::ConsIdx::ep];
    epl = L[+FieldState::ConsIdx::ep];

    Bl = L[+FieldState::ConsIdx::Bz];
    Br = R[+FieldState::ConsIdx::Bz];

    mur = R[+FieldState::ConsIdx::mu];
    mul = L[+FieldState::ConsIdx::mu];

    cr = 1/std::sqrt(mur*epr);
    cl = 1/std::sqrt(mul*epl);

    F[+FieldState::ConsIdx::Dy] = c0*((Bl*cl + Br*cr) + (Dl/epl - Dr/epr))/(cl*mul + cr*mur);
    F[+FieldState::ConsIdx::Bz] = c0*((Dl*cl + Dr*cr) + (Bl/mul - Br/mur))/(cl*epl + cr*epr);


    // div clean
    if (ch > 0) {
        Pl = L[+FieldState::ConsIdx::phi];
        Pr = R[+FieldState::ConsIdx::phi];

        fl = L[+FieldState::ConsIdx::Dx];
        fr = R[+FieldState::ConsIdx::Dx];

        F[+FieldState::ConsIdx::Dx] =   0.5*c0*(Pr + Pl)  - 0.5*ch*(fr - fl);
        F[+FieldState::ConsIdx::phi] = 0.5*ch2/c0*(fr + fl) - 0.5*ch*(Pr - Pl);
    }

    // B-wave

    Dl = L[+FieldState::ConsIdx::Dz];
    Dr = R[+FieldState::ConsIdx::Dz];

    epr = R[+FieldState::ConsIdx::ep];
    epl = L[+FieldState::ConsIdx::ep];

    Bl = L[+FieldState::ConsIdx::By];
    Br = R[+FieldState::ConsIdx::By];

    F[+FieldState::ConsIdx::Dz] = c0*(-(Bl*cl + Br*cr) + (Dl/epl - Dr/epr))/(cl*mul + cr*mur);
    F[+FieldState::ConsIdx::By] = c0*(-(Dl*cl + Dr*cr) + (Bl/mul - Br/mur))/(cl*epl + cr*epr);

    // div clean
    if (ch > 0) {
        Pl = L[+FieldState::ConsIdx::psi];
        Pr = R[+FieldState::ConsIdx::psi];

        fl = L[+FieldState::ConsIdx::Bx];
        fr = R[+FieldState::ConsIdx::Bx];

        F[+FieldState::ConsIdx::Bx] =   0.5*c0*(Pr + Pl)  - 0.5*ch*(fr - fl);
        F[+FieldState::ConsIdx::psi] = 0.5*ch2/c0*(fr + fl) - 0.5*ch*(Pr - Pl);
    }


    return;
}

bool FieldRH::valid_state(const int idx)
{

    if (MFP::get_state(idx).get_type() != State::StateType::Field) {
        return false;
    }
    return true;
}
