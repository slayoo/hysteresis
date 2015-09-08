#pragma once 

#include <blitz/array.h>

#include <libcloudph++/common/moist_air.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/maxwell-mason.hpp>
#include <libcloudph++/common/kappa_koehler.hpp>

// <TODO>
// parameters
auto p0 = 1e5 * si::pascals;
auto T0 = 3e2 * si::kelvins;
auto r0 = 2.2e-2 * si::dimensionless();
auto w  = 3e0 * si::metres_per_second;
auto N  = 1e2 * 1e6 / si::cubic_metres; // TODO: DV(rho)!
auto kpa = .8 * si::dimensionless();
auto rd3 = pow(1e-7, 3) * si::cubic_metres;
// </TODO>

struct odeset_t
{
  static constexpr double 
    t_scale = 10, // [sec] - timescale used to initialise solver
    reltol = 1e-8,
    abstol = 0;

  enum { n_eqns = 3 };
  enum class ret_t { OK, KO, adj };
  enum { ix_rv, ix_th, ix_rw };

// derived params
decltype(T0) th0;

  ret_t operator()(
    const double &t, 
    const blitz::Array<double, 1> &Y,
    blitz::Array<double, 1> dY_dt
  ) {
    namespace cmn = libcloudphxx::common;
    namespace si = boost::units::si;

    // current state with units
    auto rv = Y(ix_rv) * si::dimensionless();
    auto th = Y(ix_th) * si::kelvins;
    auto rw = Y(ix_rw) * si::metres;

    // derived quantities
    auto z = (t*si::seconds) * w;
    auto rhod = cmn::theta_std::rhod(
      cmn::hydrostatic::p(z, th0, r0, 0. * si::metres, p0),
      th0,
      r0
    );
    auto T = cmn::theta_dry::T(th, rhod);
    auto p = cmn::theta_dry::p(rhod, rv, T);
    auto p_v = rhod * rv * cmn::moist_air::R_v<double>() * T;
    auto RH = p_v / cmn::const_cp::p_vs(T);

std::cout << "t= " << t << " RH= " << RH << " rw= " << rw << " rv=" << rv << std::endl;

    // values of the derivatives
    dY_dt(ix_rw) = cmn::maxwell_mason::rdrdt(
      cmn::moist_air::D_0<double>(), // TODO: ...
      cmn::moist_air::K_0<double>(), // TODO: ...
      rhod * rv,
      T,
      p,
      RH,
      cmn::kappa_koehler::a_w(rw*rw*rw, rd3, kpa),
      cmn::kelvin::klvntrm(rw, T)
    ) / rw / (si::metres / si::second);
    
    dY_dt(ix_rv) = - 4./3. /* * pi<double>() */ * cmn::moist_air::rho_w<double>() / rhod * N 
      * 3. * rw * rw
      * (dY_dt(ix_rw) * si::metres_per_second)
      * (si::seconds);

    dY_dt(ix_th) = cmn::theta_dry::d_th_d_rv(T, th) 
      * (dY_dt(ix_rv) / si::seconds)
      / (si::kelvins / si::seconds);

//std::cerr << dY_dt << std::endl;

    return ret_t::OK;
  }

  blitz::Array<double, 1> Y0;

  // ctor
  odeset_t()
  {
    Y0.resize(n_eqns);

    // initial condition
    namespace cmn = libcloudphxx::common;

    th0 = T0 * pow(cmn::theta_std::p_1000<double>() / p0, cmn::moist_air::R_d<double>() / cmn::moist_air::c_pd<double>());
    auto RH0 = p0 * r0 / (r0 + cmn::moist_air::eps<double>()) / cmn::const_cp::p_vs(T0);

    Y0(ix_th) = cmn::theta_dry::std2dry(th0, r0) / si::kelvins;
    Y0(ix_rv) = r0;
    Y0(ix_rw) = cbrt(cmn::kappa_koehler::rw3_eq(rd3, kpa, RH0, T0) / si::cubic_metres);
  }
};
