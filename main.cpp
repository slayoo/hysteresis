#include <algorithm>
#include <limits>

#include "src/solver.hpp"
#include "src/odeset.hpp"
#include "src/invfft.hpp"                                            
#include "src/zintrp.hpp"
#include "src/hdf5io.hpp"

template <typename xi_t>
void run(hdf5io_t &hdf5io)
{
  using odeset_t = odeset_t<xi_t>;
  struct params_t : 
    zintrp_t::params_t,
    solver_t<odeset_t>::params_t,
    invfft_t::params_t,
    odeset_t::params_t
  {
    params_t(const hdf5io_t &io) :
      zintrp_t::params_t({
	.n_cycl = (int)io.getpar("n_cycl"), 
      }),
      solver_t<odeset_t>::params_t({
	.reltol = io.getpar("reltol"), 
	.abstol = io.getpar("abstol")
      }),
      invfft_t::params_t({
	.z_hlf  = io.getpar("z_hlf") * si::metres,
	.t_hlf  = io.getpar("t_hlf") * si::seconds,
	.freq   = io.getpar("freq")  * si::hertz,
	.ampl   = io.getpar("ampl")  * si::metres
      }),
      odeset_t::params_t({
	.p0     = io.getpar("p0") * si::pascals,
	.T0     = io.getpar("T0") * si::kelvins,
	.RH0     = io.getpar("RH0") * si::dimensionless(),
	.N_stp  = io.getpar("N_stp") / si::cubic_metres,
	.kpa    = io.getpar("kpa") * si::dimensionless(),
	.rd3    = pow(io.getpar("rd"), 3) * si::cubic_metres,
        .kelvin = (bool)io.getpar("kelvin"),
        .raoult = (bool)io.getpar("raoult")
      })
    {}
  } params(hdf5io);

  auto zintrp = zintrp_t(invfft_t()(params), params);
  auto odeset = odeset_t(zintrp, params);
  auto solver = solver_t<odeset_t>(odeset, params);

  while (solver.t < params.n_cycl * 2 * (params.t_hlf / si::seconds))
  {
    solver.step();

    std::array<double, hdf5io_t::n_vars> rec;
    rec[hdf5io_t::ix_t ] = solver.t;
    rec[hdf5io_t::ix_z ] = zintrp.z(solver.t);
    auto diag = odeset.diag(solver.state, solver.t);
    rec[hdf5io_t::ix_RH] = diag.RH;
    rec[hdf5io_t::ix_T] = diag.T;
    rec[hdf5io_t::ix_kelvin] = diag.kelvin;
    rec[hdf5io_t::ix_raoult] = diag.raoult;
    rec[hdf5io_t::ix_rw] = diag.rw;
    hdf5io.putrec(rec);
  }
}

int main(int argc, char **argv)
{
  if (!(argc == 2 || argc == 3))
    throw std::runtime_error("expecting one or two arguments - HDF5 filename & xi variable choice (id,p2,p3,ln)");

  hdf5io_t hdf5io(argv[1]);
  
  if (argc == 2 || std::string(argv[2]) == "id")
    run<xi_id<>>(hdf5io);
  else if (std::string(argv[2]) == "p2")
    run<xi_p2<>>(hdf5io);
  else if (std::string(argv[2]) == "p3")
    run<xi_p3<>>(hdf5io);
  else if (std::string(argv[2]) == "ln")
    run<xi_ln<>>(hdf5io);
  else
    throw std::runtime_error("unknown xi choice (second argument)");
}
