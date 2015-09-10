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
    static int ode(double t, N_Vector Y, N_Vector dY_dt, void *odeset)
    {
#if !defined(NDEBUG)
      nvec2bltz(dY_dt) = blitz::signalling_NaN(44.);
#endif
      auto ret = static_cast<odeset_t*>(odeset)->step(t, nvec2bltz(Y), nvec2bltz(dY_dt)); 
      return 
        ret == odeset_t::ret_t::OK  ? 0 : //  0: ok
        ret == odeset_t::ret_t::adj ? 1 : //  1: adjust-timestep
        -1;                               // -1: fatal
    }
  };

  void* cvode_mem;
  N_Vector Y;
  odeset_t odeset;

  public:

  double t = 0.;
  decltype(nvec2bltz(Y)) state;

  struct params_t
  {
    const double reltol, abstol;
  };

  // ctor
  solver_t(odeset_t &odeset, const params_t &params) 
    : odeset(odeset)
  {
    Y = N_VNew_Serial(odeset_t::n_eqns);
    state.reference(nvec2bltz(Y));
    odeset.init(state);

    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    if (cvode_mem == NULL)
      throw std::runtime_error("CVodeCreate() failed.");

    if (CV_SUCCESS != CVodeSetUserData(cvode_mem, &odeset))
      throw std::runtime_error("CVodeSetUserData() failed.");

    if (CV_SUCCESS != CVodeInit(cvode_mem, &(system_t::ode), t, Y))
      throw std::runtime_error("CVodeInit() failed.");

    if (CV_SUCCESS != CVodeSStolerances(cvode_mem, params.reltol, params.abstol))
      throw std::runtime_error("CVodeSStolerances() failed.");

    if (CVDLS_SUCCESS != CVDense(cvode_mem, odeset_t::n_eqns)) 
      throw std::runtime_error("CVDens() failed.");
  }

  void step()
  {
    if (CV_SUCCESS != CVode(cvode_mem, odeset_t::t_scale, Y, &t, CV_ONE_STEP)) 
      throw std::runtime_error("CVode() failed.");
  }

  // dtor
  ~solver_t() 
  {
    N_VDestroy(Y);
    CVodeFree(&cvode_mem); 
  }
};
