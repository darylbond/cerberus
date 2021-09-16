#include "MFP_hybrid_hll.H"
#include "MFP_utility.H"
#include "MFP.H"
#include "MFP_state.H"

//================================================================================

std::string HydroHybridHLL::tag = "HLLE/HLLC";
bool HydroHybridHLL::registered = GetHydroRiemannSolverFactory().Register(HydroHybridHLL::tag, HydroRiemannSolverBuilder<HydroHybridHLL>);

HydroHybridHLL::HydroHybridHLL() {}
HydroHybridHLL::HydroHybridHLL(const int i)
{
    idx = i;
    hllc = HydroHLLC(i);
    hlle = HydroHLLE(i);
}

void HydroHybridHLL::solve(Array<Real,+HydroDef::FluxIdx::NUM> &L,
                                Array<Real,+HydroDef::FluxIdx::NUM> &R,
                                Array<Real,+HydroDef::ConsIdx::NUM> &F,
                                Real* shk) const
{
    BL_PROFILE("HydroHybridHLL::solve");
    const Real eps = 1e-14;

    Array<Real,+HydroDef::ConsIdx::NUM> F_hlle, F_hllc;

    if (*shk < eps) {
        hllc.solve(L, R, F, shk);
    } else if ((1.0 - *shk) < eps) {
        hlle.solve(L, R, F, shk);
    } else {
        hllc.solve(L, R, F_hllc, shk);
        hlle.solve(L, R, F_hlle, shk);
        for (int i=0; i<+HydroDef::ConsIdx::NUM; ++i) {
            F[i] = (1.0 - *shk)*F_hllc[i] + *shk*F_hlle[i];
        }
    }
}
