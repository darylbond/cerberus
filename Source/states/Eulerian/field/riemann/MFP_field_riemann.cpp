#include "MFP_field_riemann.H"
#include <math.h>
#include <algorithm>
#include <iostream>

#include "MFP_utility.H"

FieldRiemannSolver::FieldRiemannSolver()
{
}

FieldRiemannSolver::~FieldRiemannSolver()
{
    // do nothing
}

ClassFactory<FieldRiemannSolver>& GetRiemannSolverFactory()
{
    static ClassFactory<FieldRiemannSolver> F;
    return F;
}
