#ifndef MFP_RK4_H
#define MFP_RK4_H

#include "AMReX_REAL.H"
#include "AMReX_Array.H"
#include "AMReX_Vector.H"

using amrex::Real;
using amrex::Array;
using amrex::Vector;

/*
  Function to advance set of coupled first-order o.d.e.s by single step
  using fixed step-length fourth-order Runge-Kutta scheme

     x        ... independent variable
     y        ... array of dependent variables
     h        ... fixed step-length

  Requires right-hand side routine

     void rhs_eval (Real t, Array<Real,1> y, Array<Real,1>& dydt, Real* data)

  which evaluates derivatives of y (w.r.t. t) in array dydt
*/

template<size_t N, size_t M>
int rk4_fixed (Real& t, Array<Real,N>& y, Array<Real,M>& data,
                int (*rhs_eval)(Real, Array<Real,N>, Array<Real,N>&, Array<Real,M>&),
                bool (*check)(Array<Real,N>&),
                Real h)
{
  // Array y assumed to be of extent n, where n is no. of coupled
  // equations

  // Declare local arrays
  Array<Real,N> k1, k2, k3, k4, f, dydt;

  // Zeroth intermediate step
  if((*rhs_eval) (t, y, dydt, data)) return 1;
  for (int j = 0; j < N; j++)
    {
      k1[j] = h * dydt[j];
      f[j] = y[j] + k1[j] / 2.;
    }

  // First intermediate step
  if((*rhs_eval) (t + h / 2., f, dydt, data)) return 1;
  for (int j = 0; j < N; j++)
    {
      k2[j] = h * dydt[j];
      f[j] = y[j] + k2[j] / 2.;
    }

  // Second intermediate step
  if ((*rhs_eval) (t + h / 2., f, dydt, data)) return 1;
  for (int j = 0; j < N; j++)
    {
      k3[j] = h * dydt[j];
      f[j] = y[j] + k3[j];
    }

  // Third intermediate step
  if ((*rhs_eval) (t + h, f, dydt, data)) return 1;
  for (int j = 0; j < N; j++)
    {
      k4[j] = h * dydt[j];
    }

  // Actual step
  for (int j = 0; j < N; j++)
    {
      y[j] += k1[j] / 6. + k2[j] / 3. + k3[j] / 3. + k4[j] / 6.;
    }

  if ((*check)(y)) return 1;

  t += h;

  return 0;
}

/*
  Function to advance set of coupled first-order o.d.e.s by single step
  using adaptive fourth-order Runge-Kutta scheme

     x       ... independent variable
     y       ... array of dependent variables
     h       ... step-length
     t_err   ... actual truncation error per step
     acc     ... desired truncation error per step
     S       ... step-length cannot change by more than this factor from
                  step to step
     rept    ... number of step recalculations
     maxrept ... maximum allowable number of step recalculations
     h_min   ... minimum allowable step-length
     h_max   ... maximum allowable step-length
     flag    ... controls manner in which truncation error is calculated

  Requires right-hand side routine

        void rhs_eval (Real t, Array<Real,N> y, Array<Real,N>& dydt, Real* data)

  which evaluates derivatives of y (w.r.t. t) in array dydt.

  Function advances equations by single step whilst attempting to maintain
  constant truncation error per step of acc:

    flag = 0 ... error is absolute
    flag = 1 ... error is relative
    flag = 2 ... error is mixed

  If step-length falls below h_min then routine aborts
*/

template<size_t N, size_t M>
void rk4_adaptive (Real& t, Array<Real,N>& y, Array<Real,M>& data,
                   int (*rhs_eval)(Real, Array<Real,N>, Array<Real,N>&, Array<Real,M>&),
                   bool (*check)(Array<Real,N>&),
                   Real& h, Real& t_err, Real acc,
                   Real S, int& rept, int maxrept,
                   Real h_min, Real h_max, int flag,
                   Vector<std::pair<Real,Array<Real,N>>>& history)
{

  // Declare local arrays
  Array<Real,N> y0, y1;

  // Declare repetition counter
  static int count = 0;

  // Save initial data
  Real x0 = t;
  y0 = y;

  // Take full step
  if (rk4_fixed (t, y, data, rhs_eval, check, h)) {
      // step failed, reduce time step and try again

      count++;
      t = x0;
      y = y0;
      h /= S;
      rk4_adaptive (t, y, data, rhs_eval, check, h, t_err, acc,
                    S, rept, maxrept, h_min, h_max, flag, history);
    }

  // Save data
  y1 = y;

  // Restore initial data
  t = x0;
  y = y0;

  // Take two half-steps
  if (rk4_fixed (t, y, data, rhs_eval, check, h/2.0)) {
      // step failed, reduce time step and try again

      count++;
      t = x0;
      y = y0;
      h /= S;
      rk4_adaptive (t, y, data, rhs_eval, check, h, t_err, acc,
                    S, rept, maxrept, h_min, h_max, flag, history);
    }




  if (rk4_fixed (t, y, data, rhs_eval, check, h/2.0)) {
      // step failed, reduce time step and try again

      count++;
      t = x0;
      y = y0;
      h /= S;
      rk4_adaptive (t, y, data, rhs_eval, check, h, t_err, acc,
                    S, rept, maxrept, h_min, h_max, flag, history);
    }

  // Calculate truncation error
  t_err = 0.;
  Real err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error
      for (int i = 0; i < N; i++)
        {
          err = fabs (y[i] - y1[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < N; i++)
        {
          err = fabs ((y[i] - y1[i]) / y[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else
    {
      // Use mixed truncation error
      for (int i = 0; i < N; i++)
        {
          err1 = fabs ((y[i] - y1[i]) / y[i]);
          err2 = fabs (y[i] - y1[i]);
          err = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err == 0.) t_err = 1.e-15;

  // Calculate new step-length
  Real h_est = h * pow (fabs (acc / t_err), 0.2);

  // Prevent step-length from changing by more than factor S
  if (h_est / h > S)
    h *= S;
  else if (h_est / h < 1. / S)
    h /= S;
  else
    h = h_est;

  // Prevent step-length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h / fabs(h) : h;

  // Abort if step-length falls below h_min
  if (fabs(h) < h_min)
    {
      amrex::Abort("Error - |h| < hmin\n");
    }

  // If truncation error acceptable take step
  if ((t_err <= acc) || (count >= maxrept))
    {
      rept = count;
      count = 0;
    }
  // If truncation error unacceptable repeat step
  else
    {
      count++;
      t = x0;
      y = y0;
      rk4_adaptive (t, y, data, rhs_eval, check, h, t_err, acc,
                    S, rept, maxrept, h_min, h_max, flag, history);
    }

  history.push_back(std::pair(t, y));

  return;
}

#endif // MFP_RK4_H
