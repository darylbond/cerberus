#include "MFP_hllc.H"
#include "MFP_utility.H"

//================================================================================
std::string HydroHLLC::tag = "HLLC";
bool HydroHLLC::registered = GetHydroRiemannSolverFactory().Register(HydroHLLC::tag, HydroRiemannSolverBuilder<HydroHLLC>);

HydroHLLC::HydroHLLC(){}
HydroHLLC::HydroHLLC(const int i)
{
    idx = i;
}

void HydroHLLC::solve(Vector<Real> &L,
                                Vector<Real> &R,
                                Vector<Real> &F,
                                Real* shk) const
{
    BL_PROFILE("HydroHLLC::solve");

    const int n_alpha = L.size() - +HydroDef::FluxIdx::NUM;

    // get the data out of the passed in arrays
    Real rhoL = L[+HydroDef::FluxIdx::Density];
    Real uL = L[+HydroDef::FluxIdx::Xvel];
    Real vL = L[+HydroDef::FluxIdx::Yvel];
    Real wL = L[+HydroDef::FluxIdx::Zvel];
    Real pL = L[+HydroDef::FluxIdx::Prs];
    Real gamL = L[+HydroDef::FluxIdx::Gamma];
    Real nrgL = pL/(gamL - 1) + 0.5*rhoL*(uL*uL + vL*vL + wL*wL);

    Vector<Real> trL(n_alpha);
    for (int i=+HydroDef::FluxIdx::NUM; i<L.size(); ++ i) {
        trL[i-+HydroDef::FluxIdx::NUM]= L[i]*rhoL;
    }

    Real aL = sqrt(gamL*pL/rhoL);

    // get the data out of the passed in arrays
    Real rhoR = R[+HydroDef::FluxIdx::Density];
    Real uR = R[+HydroDef::FluxIdx::Xvel];
    Real vR = R[+HydroDef::FluxIdx::Yvel];
    Real wR = R[+HydroDef::FluxIdx::Zvel];
    Real pR = R[+HydroDef::FluxIdx::Prs];
    Real gamR = R[+HydroDef::FluxIdx::Gamma];
    Real nrgR = pR/(gamR-1) + 0.5*rhoR*(uR*uR + vR*vR + wR*wR);
    Real aR = sqrt(gamR*pR/rhoR);

    Vector<Real> trR(n_alpha);
    for (int i=+HydroDef::FluxIdx::NUM; i<R.size(); ++ i) {
        trR[i-+HydroDef::FluxIdx::NUM]= R[i]*rhoR;
    }

    // Calculate wave speeds S_L, S_star and S_R

    Real rho_bar = 0.5*(rhoL + rhoR);
    Real a_bar = 0.5*(aL + aR);

    Real p_star = 0.5*(pL + pR) - 0.5*(uR - uL)*rho_bar*a_bar;
    Real S_star = 0.5*(uL + uR) - 0.5*(pR - pL)/(rho_bar*a_bar);

    Real qL;
    if (p_star <= pL) {
        qL = 1.0;
    } else {
        qL = std::sqrt(1.0 + ((gamL+1.0)/(2*gamL))*(p_star/pL - 1.0));
    }

    Real S_L = uL - aL*qL;

    Real qR;
    if (p_star <= pR) {
        qR = 1.0;
    } else {
        qR = std::sqrt(1.0 + ((gamR+1.0)/(2*gamR))*(p_star/pR - 1.0));
    }

    Real S_R = uR + aR*qR;


    if (S_L >= 0.0) {
        // flux vector L
        F[+HydroDef::ConsIdx::Density]  = rhoL*uL;
        F[+HydroDef::ConsIdx::Xmom]   = rhoL*uL*uL + pL;
        F[+HydroDef::ConsIdx::Ymom]   = rhoL*uL*vL;
        F[+HydroDef::ConsIdx::Zmom]   = rhoL*uL*wL;
        F[+HydroDef::ConsIdx::Eden]   = uL*(nrgL + pL);

        for (int i=0; i<n_alpha; ++i) {
            F[+HydroDef::ConsIdx::NUM+i]  = uL * trL[i];
        }

        return;
    } else if ((S_L <= 0.0) && (0.0 <= S_star)) {
        Array<Real, +HydroDef::ConsIdx::NUM> svLs, fvL, svL;
        // flux vector L
        fvL[+HydroDef::ConsIdx::Density]  = rhoL*uL;
        fvL[+HydroDef::ConsIdx::Xmom]   = rhoL*uL*uL + pL;
        fvL[+HydroDef::ConsIdx::Ymom]   = rhoL*uL*vL;
        fvL[+HydroDef::ConsIdx::Zmom]   = rhoL*uL*wL;
        fvL[+HydroDef::ConsIdx::Eden]   = uL*(nrgL + pL);

        for (int i=0; i<n_alpha; ++i) {
            fvL[+HydroDef::ConsIdx::NUM+i]  = uL * trL[i];
        }

        // state vector L
        svL[+HydroDef::ConsIdx::Density]  = rhoL;
        svL[+HydroDef::ConsIdx::Xmom]   = rhoL*uL;
        svL[+HydroDef::ConsIdx::Ymom]   = rhoL*vL;
        svL[+HydroDef::ConsIdx::Zmom]   = rhoL*wL;
        svL[+HydroDef::ConsIdx::Eden]   = nrgL;

        for (int i=0; i<n_alpha; ++i) {
            svL[+HydroDef::ConsIdx::NUM+i]  = trL[i];
        }

        Real coeff = rhoL*((S_L - uL)/(S_L - S_star));

        svLs[+HydroDef::ConsIdx::Density] = coeff;
        svLs[+HydroDef::ConsIdx::Xmom] = coeff*S_star;
        svLs[+HydroDef::ConsIdx::Ymom] = coeff*vL;
        svLs[+HydroDef::ConsIdx::Zmom] = coeff*wL;
        svLs[+HydroDef::ConsIdx::Eden] = coeff*(nrgL/rhoL + (S_star - uL)*(S_star + pL/(rhoL*(S_L - uL))));

        for (int i=0; i<n_alpha; ++i) {
            svLs[+HydroDef::ConsIdx::NUM+i]  = trL[i]*((S_L - uL)/(S_L - S_star));
        }

        for (int i=0; i<+HydroDef::ConsIdx::NUM+n_alpha; ++i) {
            F[i] = fvL[i] + S_L*(svLs[i] - svL[i]);
        }

    } else if ((S_star <= 0.0) && (0.0 <= S_R)) {
        Array<Real, +HydroDef::ConsIdx::NUM> svRs, fvR, svR;
        // flux vector R
        fvR[+HydroDef::ConsIdx::Density]  = rhoR*uR;
        fvR[+HydroDef::ConsIdx::Xmom]   = rhoR*uR*uR + pR;
        fvR[+HydroDef::ConsIdx::Ymom]   = rhoR*uR*vR;
        fvR[+HydroDef::ConsIdx::Zmom]   = rhoR*uR*wR;
        fvR[+HydroDef::ConsIdx::Eden]   = uR*(nrgR + pR);

        for (int i=0; i<n_alpha; ++i) {
            fvR[+HydroDef::ConsIdx::NUM+i]  = uR * trR[i];
        }

        // state vector R
        svR[+HydroDef::ConsIdx::Density]  = rhoR;
        svR[+HydroDef::ConsIdx::Xmom]   = rhoR*uR;
        svR[+HydroDef::ConsIdx::Ymom]   = rhoR*vR;
        svR[+HydroDef::ConsIdx::Zmom]   = rhoR*wR;
        svR[+HydroDef::ConsIdx::Eden]   = nrgR;

        for (int i=0; i<n_alpha; ++i) {
            svR[+HydroDef::ConsIdx::NUM+i]  = trR[i];
        }

        Real coeff = rhoR*((S_R - uR)/(S_R - S_star));

        svRs[+HydroDef::ConsIdx::Density] = coeff;
        svRs[+HydroDef::ConsIdx::Xmom] = coeff*S_star;
        svRs[+HydroDef::ConsIdx::Ymom] = coeff*vR;
        svRs[+HydroDef::ConsIdx::Zmom] = coeff*wR;
        svRs[+HydroDef::ConsIdx::Eden] = coeff*(nrgR/rhoR + (S_star - uR)*(S_star + pR/(rhoR*(S_R - uR))));

        for (int i=0; i<n_alpha; ++i) {
            svR[+HydroDef::ConsIdx::NUM+i]  = trR[i]*((S_R - uR)/(S_R - S_star));
        }

        for (int i=0; i<+HydroDef::ConsIdx::NUM+n_alpha; ++i) {
            F[i] = fvR[i] + S_R*(svRs[i] - svR[i]);
        }

    } else {
        // flux vector R
        F[+HydroDef::ConsIdx::Density]  = rhoR*uR;
        F[+HydroDef::ConsIdx::Xmom]   = rhoR*uR*uR + pR;
        F[+HydroDef::ConsIdx::Ymom]   = rhoR*uR*vR;
        F[+HydroDef::ConsIdx::Zmom]   = rhoR*uR*wR;
        F[+HydroDef::ConsIdx::Eden]   = uR*(nrgR + pR);

        for (int i=0; i<n_alpha; ++i) {
            F[+HydroDef::ConsIdx::NUM+i]  = trR[i]*uR;
        }
    }
}
