#pragma once

#include <gsl/gsl_spline.h>

#include "gslerr.hpp"

#include <boost/units/systems/si.hpp>
using namespace boost::units;

class zintrp_t
{
  gsl_interp_accel *accl;
  gsl_interp *intp;

  public:

  std::array<blitz::Array<double, 1>,2> data;

  zintrp_t(std::array<blitz::Array<double, 1>,2> &&data)
    : data(data)
  {
    gslerr_init();

    if (data[0].size() != data[1].size())
      throw std::runtime_error("arrays not of same size!");

    accl = gsl_interp_accel_alloc();
    intp = gsl_interp_alloc(gsl_interp_cspline_periodic, data[0].size());
    gsl_interp_init(intp, data[0].data(), data[1].data(), data[0].size());
  }

  quantity<si::length> z(const quantity<si::time> &t)
  {
    double z;
    auto status = gsl_interp_eval_e(intp, data[0].data(), data[1].data(), t / si::seconds, accl, &z);
    assert(status == GSL_SUCCESS);
    return z * si::metres;
  }

  ~zintrp_t()
  {
    gsl_interp_free(intp);
    gsl_interp_accel_free(accl);
  }
};
