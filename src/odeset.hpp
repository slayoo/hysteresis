#pragma once 

#include <blitz/array.h>

#include <libcloudph++/common/earth.hpp>
#include <libcloudph++/common/moist_air.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/maxwell-mason.hpp>
#include <libcloudph++/common/kappa_koehler.hpp>
#include <libcloudph++/common/mean_free_path.hpp>
#include <libcloudph++/common/transition_regime.hpp>
using namespace libcloudphxx::common;

#include "zintrp.hpp"

class odeset_t
{
  public:

  struct params_t
  {
    boost::units::quantity<si::pressure> p0;
    boost::units::quantity<si::temperature> T0;
    boost::units::quantity<si::dimensionless> r0;
    decltype(1./si::cubic_metres) N_stp;
    boost::units::quantity<si::dimensionless> kpa;
    boost::units::quantity<si::volume> rd3;
  };

  private:

  // set in ctor, used in init() and step()
  quantity<si::temperature> th0;
  quantity<si::mass_density> rhod0;

  const zintrp_t &zintrp;
  const params_t &params;

  public:

  static constexpr double t_scale = 10; // [sec] - timescale used to initialise solver

  enum { n_eqns = 3 };
  enum class ret_t { OK, KO, adj };
  enum { ix_rv, ix_th, ix_rw };

  ret_t step(
    const double &t, 
    const blitz::Array<double, 1> &Y,
    blitz::Array<double, 1> dY_dt
  ) {
    namespace si = boost::units::si;

    // current state with units
    auto rv = Y(ix_rv) * si::dimensionless();
    auto th = Y(ix_th) * si::kelvins;
    auto rw = Y(ix_rw) * si::metres;

    // derived quantities
    auto rhod = _rhod(t * si::seconds);
    auto T = theta_dry::T(th, rhod);
    auto p = theta_dry::p(rhod, rv, T);
    auto p_v = rhod * rv * moist_air::R_v<double>() * T;

    // values of the derivatives
    dY_dt(ix_rw) = maxwell_mason::rdrdt(
      moist_air::D_0<double>() * transition_regime::beta(mean_free_path::lambda_D(T)    / rw), 
      moist_air::K_0<double>() * transition_regime::beta(mean_free_path::lambda_K(T, p) / rw),
      rhod * rv, T, p,
      _RH(T, rv, rhod),
      kappa_koehler::a_w(rw*rw*rw, params.rd3, params.kpa),
      kelvin::klvntrm(rw, T)
    ) / rw / (si::metres / si::second);
    
    dY_dt(ix_rv) = - 4./3. * pi<double>() * moist_air::rho_w<double>() 
      * params.N_stp * (rhod0 / earth::rho_stp<double>()) / rhod  //
      * 3. * rw * rw  // 
      * (dY_dt(ix_rw) * si::metres_per_second)
      * (si::seconds);

    dY_dt(ix_th) = theta_dry::d_th_d_rv(T, th) 
      * (dY_dt(ix_rv) / si::seconds)
      / (si::kelvins / si::seconds);

    return ret_t::OK;
  }

  // initial condition
  void init(blitz::Array<double, 1> Y)
  {
    auto RH0 = params.p0 * params.r0 / (params.r0 + moist_air::eps<double>()) / const_cp::p_vs(params.T0);

    Y(ix_th) = theta_dry::std2dry(th0, params.r0) / si::kelvins;
    Y(ix_rv) = params.r0;
    Y(ix_rw) = cbrt(kappa_koehler::rw3_eq(params.rd3, params.kpa, RH0, params.T0) / si::cubic_metres);
  }

  // ctor
  odeset_t(const zintrp_t &zintrp, const params_t &params)
    : zintrp(zintrp), params(params)
  {
    th0 = params.T0 * pow(theta_std::p_1000<double>() / params.p0, moist_air::R_d<double>() / moist_air::c_pd<double>());
    rhod0 = theta_std::rhod(params.p0, th0, params.r0);
  }

  // for diagnostics
  private:

  quantity<si::mass_density> _rhod(const quantity<si::time> &t) const
  {
    return theta_std::rhod(
      hydrostatic::p(zintrp.z(t / si::seconds) * si::metres, th0, params.r0, 0. * si::metres, params.p0),
      th0,
      params.r0
    );
  }

  quantity<si::dimensionless> _RH(
    const quantity<si::temperature> &T,
    const quantity<si::dimensionless> &rv,
    const quantity<si::mass_density> &rhod
  ) const {
    auto p_v = rhod * rv * moist_air::R_v<double>() * T;
    return p_v / const_cp::p_vs(T);
  }

  public:

  struct diag_t
  {
    double T, RH, kelvin, raoult;
  };

  diag_t diag(
    const blitz::Array<double, 1> &Y,
    const double &t
  ) const {
    auto rhod = _rhod(t * si::seconds);
    diag_t ret;
    ret.T = theta_dry::T(Y(ix_th) * si::kelvins, rhod) / si::kelvins;
    ret.kelvin = kelvin::klvntrm(Y(ix_rw) * si::metres, ret.T * si::kelvins); // TODO: Jacobian
    ret.raoult = kappa_koehler::a_w(pow(Y(ix_rw),3) * si::cubic_metres, params.rd3, params.kpa); // TODO: Jacobian
    ret.RH = _RH(ret.T * si::kelvins, Y(ix_rv) * si::dimensionless(), rhod);
    return ret;
  }
};
