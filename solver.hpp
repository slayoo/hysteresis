#pragma once 

#include <stdexcept>
#include <iostream>

#include <blitz/array.h>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

template <typename odeset_t>
class solver_t
{
  double t;

  static auto nvec2bltz(const N_Vector &nvec)
  {
    return blitz::Array<double, 1>(
      N_VGetArrayPointer(nvec),
      blitz::TinyVector<int, 1>({odeset_t::n_eqns}), 
      blitz::neverDeleteData
    );
  }

  struct system_t 
  {
    static int ode(double t, N_Vector Y, N_Vector dY_dt, void *params)
    {
      auto ret = static_cast<params_t*>(params)->odeset(t, nvec2bltz(Y), nvec2bltz(dY_dt)); 
      return 
        ret == odeset_t::ret_t::OK  ? 0 : //  0: ok
        ret == odeset_t::ret_t::adj ? 1 : //  1: adjust-timestep
        -1;                               // -1: fatal
    }
  };

  void* cvode_mem;
  N_Vector Y;

  struct params_t {
    odeset_t odeset;
  } params;

  public:

  solver_t(odeset_t &&odeset)
    : params({odeset})
  {
    Y = N_VNew_Serial(odeset_t::n_eqns);
    nvec2bltz(Y) = odeset.Y0;

    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    if (cvode_mem == NULL)
      throw std::runtime_error("CVodeCreate() failed.");

    if (CV_SUCCESS != CVodeSetUserData(cvode_mem, &params))
      throw std::runtime_error("CVodeSetUserData() failed.");

    if (CV_SUCCESS != CVodeInit(cvode_mem, &(system_t::ode), /* t0 = */ 0., Y))
      throw std::runtime_error("CVodeInit() failed.");

    if (CV_SUCCESS != CVodeSStolerances(cvode_mem, odeset_t::reltol, odeset_t::abstol))
      throw std::runtime_error("CVodeSStolerances() failed.");

    if (CVDLS_SUCCESS != CVDense(cvode_mem, odeset_t::n_eqns)) 
      throw std::runtime_error("CVDens() failed.");
  }

  void step()
  {
    if (CV_SUCCESS != CVode(cvode_mem, odeset_t::t_scale, Y, &t, CV_ONE_STEP)) 
      throw std::runtime_error("CVode() failed.");
  }

  ~solver_t() 
  {
    N_VDestroy(Y);
    CVodeFree(&cvode_mem); 
  }
};
