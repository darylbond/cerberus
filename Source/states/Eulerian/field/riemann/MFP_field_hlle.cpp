#include "MFP_field_hlle.H"
#include "MFP_utility.H"
#include "MFP.H"
#include "MFP_field.H"

//================================================================================

std::string FieldHLLE::tag = "HLLE";
bool FieldHLLE::registered = GetRiemannSolverFactory().Register(FieldHLLE::tag, RiemannSolverBuilder<FieldHLLE>);

FieldHLLE::FieldHLLE(){}
FieldHLLE::FieldHLLE(const int i)
{
    idx = i;

    FieldState& istate = static_cast<FieldState&>(MFP::get_state(i));

    Real c0 = MFP::lightspeed;
    Real c2 = c0*c0;
    Real ch = istate.div_speed;
    Real ch2 = ch*ch;
    Real cc = ch2/c2;
}

void FieldHLLE::solve(Array<Real,+FieldState::ConsIdx::NUM> &L,
                      Array<Real,+FieldState::ConsIdx::NUM> &R,
                      Array<Real,+FieldState::ConsIdx::NUM> &F) const
{
    BL_PROFILE("FieldHLLE::solve");

    Array<Real, +FieldState::ConsIdx::NUM> FL, FR;

    FL[+FieldState::ConsIdx::Bx]   =   L[+FieldState::ConsIdx::psi];
    FL[+FieldState::ConsIdx::By]   = - L[+FieldState::ConsIdx::Dz]/L[+FieldState::ConsIdx::ep];
    FL[+FieldState::ConsIdx::Bz]   =   L[+FieldState::ConsIdx::Dy]/L[+FieldState::ConsIdx::ep];
    FL[+FieldState::ConsIdx::Dx]   =   L[+FieldState::ConsIdx::phi];
    FL[+FieldState::ConsIdx::Dy]   =   L[+FieldState::ConsIdx::Bz]/L[+FieldState::ConsIdx::mu];
    FL[+FieldState::ConsIdx::Dz]   = - L[+FieldState::ConsIdx::By]/L[+FieldState::ConsIdx::mu];
    FL[+FieldState::ConsIdx::psi] =   L[+FieldState::ConsIdx::Bx]*cc;
    FL[+FieldState::ConsIdx::phi] =   L[+FieldState::ConsIdx::Dx]*cc;

    FR[+FieldState::ConsIdx::Bx]   =   R[+FieldState::ConsIdx:: psi];
    FR[+FieldState::ConsIdx::By]   = - R[+FieldState::ConsIdx:: Dz]/R[+FieldState::ConsIdx::ep];
    FR[+FieldState::ConsIdx::Bz]   =   R[+FieldState::ConsIdx:: Dy]/R[+FieldState::ConsIdx::ep];
    FR[+FieldState::ConsIdx::Dx]   =   R[+FieldState::ConsIdx:: phi];
    FR[+FieldState::ConsIdx::Dy]   =   R[+FieldState::ConsIdx:: Bz]/R[+FieldState::ConsIdx::mu];
    FR[+FieldState::ConsIdx::Dz]   = - R[+FieldState::ConsIdx:: By]/R[+FieldState::ConsIdx::mu];
    FR[+FieldState::ConsIdx::psi] =   R[+FieldState::ConsIdx::Bx]*cc;
    FR[+FieldState::ConsIdx::phi] =   R[+FieldState::ConsIdx::Dx]*cc;

    Real speed = c0; // UPDATE THIS TO ACCOMMODATE DIV SPEED

    for (int n=0; n<+FieldState::ConsIdx::NUM; ++n) {
        F[n] = 0.5*(FL[n]*c0 + FR[n]*c0 + (L[n] - R[n])*speed);
    }

}

bool FieldHLLE::valid_state(const int idx)
{

    if (MFP::get_state(idx).get_type() != State::StateType::Field) {
        return false;
    }
    return true;
}
