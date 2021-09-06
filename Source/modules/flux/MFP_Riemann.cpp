#include "MFP_Riemann.H"
#include <math.h>
#include <algorithm>
#include <iostream>

#include "MFP_utility.H"

RiemannSolver::RiemannSolver()
{
}

RiemannSolver::~RiemannSolver()
{
    // do nothing
}

void RiemannSolver::solve(Vector<Real> &L,
                               Vector<Real> &R,
                               Vector<Real> &F,
                               Real* shk)
{
    // do nothing
}

ClassFactory<RiemannSolver>& GetRiemannSolverFactory()
{
    static ClassFactory<RiemannSolver> F;
    return F;
}
