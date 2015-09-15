#include <algorithm>
#include <limits>

#include "src/solver.hpp"
#include "src/odeset.hpp"
#include "src/invfft.hpp"                                            
#include "src/zintrp.hpp"
#include "src/hdf5io.hpp"

int main(int argc, char **argv)
{
  if (argc != 2)
    throw std::runtime_error("expecting one argument - HDF5 filename");

  hdf5io_t hdf5io(argv[1]);

  struct params_t : 
    solver_t<odeset_t>::params_t,
    invfft_t::params_t,
    odeset_t::params_t
  {
    params_t(const hdf5io_t &io) :
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
	.r0     = io.getpar("r0") * si::dimensionless(),
	.N_stp  = io.getpar("N_stp") / si::cubic_metres,
	.kpa    = io.getpar("kpa") * si::dimensionless(),
	.rd3    = pow(io.getpar("rd"), 3) * si::cubic_metres
      })
    {}
  } params(hdf5io);

  auto zintrp = zintrp_t(invfft_t()(params));
  auto odeset = odeset_t(zintrp, params);
  auto solver = solver_t<odeset_t>(odeset, params);

  while (solver.t < 2. * params.t_hlf / si::seconds)
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
    rec[hdf5io_t::ix_rw] = solver.state(odeset_t::ix_rw);
    hdf5io.putrec(rec);
  }
}
