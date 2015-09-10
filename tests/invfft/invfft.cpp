#include <set>
#include <boost/units/io.hpp>

#include "invfft.hpp"

// intentionally after invfft.hp so Blitz++ support is loaded
#include <gnuplot-iostream.h>

int main()
{
  Gnuplot gp;

  gp << "set term svg dynamic enhanced mouse standalone fsize 18\n";
  gp << "set output 'invfft.svg'\n";
  gp << "set grid\n";
  gp << "set key invert samplen .5\n";
  gp << "set xlabel 'time [s]'\n";
  gp << "set ylabel 'displacement [m]'\n";
  

  using arg_t = std::tuple<
    quantity<si::length>,    // z_hlf
    quantity<si::time>,      // t_hlf
    quantity<si::frequency>, // freq
    quantity<si::length>,    // ampl
    std::string              // gnuplot lt
  >;

  std::set<arg_t> cases{
    arg_t{100. * si::metres, 100. * si::seconds,  0.   * si::hertz, 0.  * si::metres, "lc rgb 'black'" },
    arg_t{100. * si::metres, 100. * si::seconds,   .05 * si::hertz, 5.  * si::metres, "lc rgb 'orange'" },
    arg_t{ 80. * si::metres, 150. * si::seconds,  0.   * si::hertz, 0.  * si::metres, "lc rgb 'blue'" },
    arg_t{ 80. * si::metres, 150. * si::seconds,   .25 * si::hertz, 1. * si::metres,  "lc rgb 'red'" },
  };

  gp << "plot 1./0 notitle "; // for the commas below :)
  for (auto &arg : cases)
  {
    std::ostringstream tmp;
    tmp << std::setprecision(3);
    tmp << "f=" << double(std::get<2>(arg) / si::hertz) << " Hz";
    tmp << "   ";
    tmp << "A=" << double(std::get<3>(arg) / si::metres) << " m";
    gp << ", '-' with lines " << std::get<4>(arg) << " title '" << tmp.str() << "'";
  }
  gp << "\n";
  for (auto &arg : cases)
    gp.send(invfft(
      std::get<0>(arg),
      std::get<1>(arg),
      std::get<2>(arg),
      std::get<3>(arg)
    ));
}
