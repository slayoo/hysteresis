#include "invfft.hpp"
#include "zintrp.hpp"

// intentionally after above so Blitz++ support is loaded
#include <gnuplot-iostream.h>

int main()
{
  namespace si = boost::units::si;

  const double xmax = .05;
  const int npts = 64;

  Gnuplot gp;

  gp << "set term svg dynamic enhanced mouse standalone fsize 18\n";
  gp << "set output 'zintrp.svg'\n";
  gp << "set grid\n";
  gp << "set key Left reverse samplen .5\n";
  gp << "set xlabel 'time [s]'\n";
  gp << "set ylabel 'displacement [m]'\n";

  gp << "set xrange [0:" << xmax << "]\n";
  gp << "set yrange [0:.01]\n";
  
  zintrp_t zintrp(invfft_t()(invfft_t::params_t({
    .z_hlf = 1.  * si::metres,  // z_hlf 
    .t_hlf = 1.  * si::seconds, // t_hlf
    .freq  = 1.  * si::hertz,   // freq
    .ampl  = .1 * si::metres   // ampl
  })));

  std::remove_const<decltype(zintrp.data)>::type dense;
  for (auto &arr : dense) arr.resize(npts);

  for (int i=0; i < npts; ++i)
  {
    dense[0](i) = double(i)/(npts-1) * xmax;
    dense[1](i) = zintrp.z(dense[0](i));
  }

  gp << "plot "
    << "'-' with linesp t 'FFT result',"
    << "'-' with linesp ps .5 pt 6 title 'spline-interpolated values'\n";
  gp.send(zintrp.data);
  gp.send(dense);
}
