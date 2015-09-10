#pragma once

#include <gsl/gsl_spline.h>

#include "gslerr.hpp"

class zintrp_t
{
  gsl_interp_accel *accl;
  gsl_interp *intp;

  public:

  const std::array<blitz::Array<double, 1>,2> data;

  // ctor
  zintrp_t(std::array<blitz::Array<double, 1>,2> &&data)
    : data(std::move(data))
  {
    gslerr_init();

    if (data[0].size() != data[1].size())
      throw std::runtime_error("arrays not of same size!");

    accl = gsl_interp_accel_alloc();
    intp = gsl_interp_alloc(gsl_interp_cspline_periodic, data[0].size());
    int status;
    status = gsl_interp_init(intp, data[0].data(), data[1].data(), data[0].size());
    assert(status == GSL_SUCCESS);
  }

  double z(
    const double &t
  ) const
  {
    if (t > data[0](data[0].size()-1)) return 0.;

    double z;
    auto status = gsl_interp_eval_e(
      intp, 
      data[0].data(), 
      data[1].data(), 
      t,
      accl, 
      &z
    );
    assert(status == GSL_SUCCESS);
    return z;
  }

  // dtor
  ~zintrp_t()
  {
    gsl_interp_free(intp);
    gsl_interp_accel_free(accl);
  }

/*
  // disabling copy ctor
  zintrp_t(const zintrp_t&) = delete;

  // enabling move ctor
  zintrp_t(zintrp_t&&) = default;
*/
};
