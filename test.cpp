#include <algorithm>
#include <limits>

#include "src/solver.hpp"
#include "src/odeset.hpp"
#include "src/invfft.hpp"                                            
#include "src/zintrp.hpp"
#include "src/hdf5io.hpp"

int main()
{
  struct params_t : 
    solver_t<odeset_t>::params_t,
    invfft_t::params_t,
    odeset_t::params_t
  {
    params_t() :
      solver_t<odeset_t>::params_t({
	.reltol = 1e-5, 
	.abstol = 0.
      }),
      invfft_t::params_t({
	.z_hlf  = 200. * si::metres,
	.t_hlf  = 200. * si::seconds,
	.freq   = .5 * si::hertz,
	.ampl   = 1. * si::metres
      }),
      odeset_t::params_t({
	.p0     = 1e5 * si::pascals,
	.T0     = 3e2 * si::kelvins,
	.r0     = 2.2e-2 * si::dimensionless(),
	.N_stp  = 1e2 * 1e6 / si::cubic_metres,
	.kpa    = .5 * si::dimensionless(),
	.rd3    = pow(1e-7, 3) * si::cubic_metres
      })
    {}
  } params;

  auto zintrp = zintrp_t(invfft_t()(params));
  auto odeset = odeset_t(zintrp, params);
  auto solver = solver_t<odeset_t>(odeset, params);

  auto t_last = solver.t, dtmin = std::numeric_limits<double>::infinity();
  while (solver.t < 2. * params.t_hlf / si::seconds)
  {
    solver.step();
    auto RH = odeset.RH(solver.state, solver.t);
std::cout 
  << solver.t
  << " " 
  << zintrp.z(solver.t)
  << " " 
  << RH
  << std::endl;
    if (RH > 1) dtmin = std::min(dtmin, solver.t - t_last);
    t_last = solver.t;
  }
std::cerr << dtmin <<  std::endl;
}
