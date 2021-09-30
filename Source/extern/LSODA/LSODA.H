/***
 *        Created:  2018-08-14

 *         Author:  Dilawar Singh <dilawars@ncbs.res.in>
 *   Organization:  NCBS Bangalore
 *        License:  MIT License
 */

#ifndef LSODE_H
#define LSODE_H

#include <array>
#include <cmath>
#include <memory>

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

using namespace std;
using amrex::Real;
using amrex::Vector;


/* --------------------------------------------------------------------------*/
/**
 * @Synopsis  Type definition of LSODA ode system.
 *
 * @Param time, Real
 * @Param y, array of Real.
 * @Param dydt, array of Real
 * @Param data, void*
 *
 * @Returns void
 */
/* ----------------------------------------------------------------------------*/

typedef int (*LSODA_ODE_SYSTEM_TYPE)(double t, Vector<Real> &y, Vector<Real> &dydt, const int offset);

class LSODA
{

public:
    LSODA();
    ~LSODA();

    size_t idamax1(const Vector<Real> &dx, const size_t n_, const size_t offset);

    void dscal1(const Real da, Vector<Real> &dx, const size_t n_,
                const size_t offset);

    Real ddot1(const Vector<Real> &a, const Vector<Real> &b, const size_t n_,
                 const size_t offsetA, const size_t offsetB);

    void daxpy1(const Real da, const Vector<Real> &dx, Vector<Real> &dy,
                const size_t n_, const size_t offsetX,
                const size_t offsetY);

    void dgesl(const Vector<Vector<Real>> &a, const size_t n_, Vector<int> &ipvt_,
               Vector<Real> &b, const size_t job);

    void dgefa(Vector<Vector<Real>> &a, const size_t n_, Vector<int> &ipvt_,
               size_t *const info);

    void prja(const size_t neq, Vector<Real> &y, LSODA_ODE_SYSTEM_TYPE f);

    void lsoda(LSODA_ODE_SYSTEM_TYPE f, const size_t neq, Vector<Real> &y,
               Real *t, Real tout, int itask, int *istate, int iopt, int jt,
               array<int, 7> &iworks, array<Real, 4> &rworks);

    void correction(const size_t neq, Vector<Real> &y, LSODA_ODE_SYSTEM_TYPE f,
                    size_t *corflag, Real pnorm, Real *del, Real *delp,
                    Real *told, size_t *ncf, Real *rh, size_t *m);

    void stoda(const size_t neq, Vector<Real> &y, LSODA_ODE_SYSTEM_TYPE f);

    // We call this function in VoxelPools::
    void lsoda_update(LSODA_ODE_SYSTEM_TYPE f, const size_t neq,
                      Vector<Real> &y, Vector<Real> &yout, Real *t,
                      const Real tout, int *istate,
                      Real rtol = 1e-6, Real atol = 1e-6, // Tolerance
                      int max_step=500);

    void terminate(int *istate);
    void terminate2(Vector<Real> &y, Real *t);
    void successreturn(Vector<Real> &y, Real *t, int itask, int ihit,
                       Real tcrit, int *istate);
    void _freeVectors(void);
    void ewset(const Vector<Real> &ycur);
    void resetcoeff(void);
    void solsy(Vector<Real> &y);
    void endstoda(void);
    void orderswitch(Real *rhup, Real dsm, Real *pdh, Real *rh,
                     size_t *orderflag);
    void intdy(Real t, int k, Vector<Real> &dky, int *iflag);
    void corfailure(Real *told, Real *rh, size_t *ncf, size_t *corflag);
    void methodswitch(Real dsm, Real pnorm, Real *pdh, Real *rh);
    void cfode(int meth);
    void scaleh(Real *rh, Real *pdh);
    Real fnorm(int n_, const Vector<Vector<Real>> &a, const Vector<Real> &w);
    Real vmnorm(const size_t n_, const Vector<Real> &v,
                  const Vector<Real> &w);

    static bool abs_compare(Real a, Real b);

private:
    size_t ml, mu, imxer;
    Real sqrteta;

    // NOTE: initialize in default constructor. Older compiler e.g. 4.8.4 would
    // produce error if these are initialized here. With newer compiler,
    // initialization can be done here.
    array<size_t, 3> mord;
    array<Real, 13> sm1;

    array<Real, 14> el;  // = {0};
    array<Real, 13> cm1; // = {0};
    array<Real, 6> cm2;  // = {0};

    array<array<Real, 14>, 13> elco;
    array<array<Real, 4>, 13> tesco;

    size_t illin, init, ierpj, iersl, jcur, l, miter, maxord, maxcor, msbp, mxncf;

    int kflag, jstart;

    size_t ixpr = 0, jtyp, mused, mxordn, mxords = 12;
    size_t meth_;

    size_t n, nq, nst, nfe, nje, nqu;
    size_t mxstep, mxhnil;
    size_t nslast, nhnil, ntrep, nyh;

    Real ccmax, el0, h_ = .0;
    Real hmin, hmxi, hu, rc, tn_ = 0.0;
    Real tsw, pdnorm;
    Real conit, crate, hold, rmax;

    size_t ialth, ipup, lmax;
    size_t nslp;
    Real pdest, pdlast, ratio;
    int icount, irflag;

    Vector<Real> ewt;
    Vector<Real> savf;
    Vector<Real> acor;
    Vector<Vector<Real>> yh_;
    Vector<Vector<Real>> wm_;

    Vector<int> ipvt;

private:
    int itol_ = 2;
    Vector<Real> rtol_;
    Vector<Real> atol_;

public:
    void *param = nullptr;
};

#endif /* end of include guard: LSODE_H */

