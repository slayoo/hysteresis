#include <algorithm>
#include <limits>

#include "src/solver.hpp"
#include "src/odeset.hpp"
#include "src/invfft.hpp"                                            
#include "src/zintrp.hpp"

struct params_t : 
  solver_t<odeset_t>::params_t,
  invfft_t::params_t
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
    })
  {}
} params;

int main()
{
  auto zintrp = zintrp_t(invfft_t()(params));
  auto odeset = odeset_t(zintrp);
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
