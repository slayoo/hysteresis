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

// <TODO>
// parameters
auto p0 = 1e5 * si::pascals;
auto T0 = 3e2 * si::kelvins;
auto r0 = 2.2e-2 * si::dimensionless();
auto w  = 3e0 * si::metres_per_second;
auto N_stp  = 1e2 * 1e6 / si::cubic_metres;
auto kpa = .5 * si::dimensionless();
auto rd3 = pow(1e-7, 3) * si::cubic_metres;
// </TODO>

class odeset_t
{
  // set in ctor, used in init() and step()
  quantity<si::temperature> th0;
  quantity<si::mass_density> rhod0;

  public:

  static constexpr double 
    t_scale = 10, // [sec] - timescale used to initialise solver
    reltol = 1e-5,
    abstol = 0;

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
      kappa_koehler::a_w(rw*rw*rw, rd3, kpa),
      kelvin::klvntrm(rw, T)
    ) / rw / (si::metres / si::second);
    
    dY_dt(ix_rv) = - 4./3. * pi<double>() * moist_air::rho_w<double>() 
      * N_stp * (rhod0 / earth::rho_stp<double>()) / rhod  //
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

    auto RH0 = p0 * r0 / (r0 + moist_air::eps<double>()) / const_cp::p_vs(T0);

    Y(ix_th) = theta_dry::std2dry(th0, r0) / si::kelvins;
    Y(ix_rv) = r0;
    Y(ix_rw) = cbrt(kappa_koehler::rw3_eq(rd3, kpa, RH0, T0) / si::cubic_metres);
  }

  // ctor
  odeset_t()
  {
    th0 = T0 * pow(theta_std::p_1000<double>() / p0, moist_air::R_d<double>() / moist_air::c_pd<double>());
    rhod0 = theta_std::rhod(p0, th0, r0);
  }

  // for diagnostics
  private:

  quantity<si::mass_density> _rhod(const quantity<si::time> &t)
  {
    auto z = t * w;
    return theta_std::rhod(
      hydrostatic::p(z, th0, r0, 0. * si::metres, p0),
      th0,
      r0
    );
  }

  quantity<si::dimensionless> _RH(
    const quantity<si::temperature> &T,
    const quantity<si::dimensionless> &rv,
    const quantity<si::mass_density> &rhod
  ) {
    auto p_v = rhod * rv * moist_air::R_v<double>() * T;
    return p_v / const_cp::p_vs(T);
  }

  public:

  double RH(
    const blitz::Array<double, 1> &Y,
    const double &t
  ) {
    auto rhod = _rhod(t * si::seconds);
    return _RH(
      theta_dry::T(Y(ix_th) * si::kelvins, rhod),
      Y(ix_rv) * si::dimensionless(),
      rhod
    );
  }
};
