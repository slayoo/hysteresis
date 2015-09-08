#pragma once 

#include <blitz/array.h>

struct odeset_t
{
  static constexpr double 
    t_scale = 10, // [sec] - timescale used to initialise solver
    reltol = 1e-8,
    abstol = 0;

  enum { n_eqns = 4 };
  enum class ret_t { OK, KO, adj };

  ret_t operator()(
    const double &t, 
    const blitz::Array<double, 1> &Y,
    blitz::Array<double, 1> dY_dt
  ) {
    // values of the derivatives
    dY_dt = 1;
    return ret_t::OK;
  }

  blitz::Array<double, 1> Y0;

  // ctor
  odeset_t()
  {
    Y0.resize(n_eqns);
    // initial condition
    Y0 = 1,2,3,4;
  }
};
