#include <array>
#include <gsl/gsl_fft_halfcomplex.h>
#include <blitz/array.h>
#include <boost/units/systems/si.hpp>

#include "gslerr.hpp"

using namespace boost::units;

std::array<blitz::Array<double, 1>, 2> invfft(
  const quantity<si::length> &z_hlf, //= 100. * si::metres;
  const quantity<si::time> &t_hlf, //= 100. * si::seconds;
  const quantity<si::frequency> &freq, //= 0. * si::hertz;
  const quantity<si::length> &ampl //= 1. * si::metres;
)
{
  gslerr_init();

  // at least 64 points per cycle
  const int n = pow(2, ceil(log2(64. * (
    freq == 0. * si::hertz 
      ? 1.
      : (freq * (2. * t_hlf)).value()
  ))));

  std::array<blitz::Array<double, 1>, 2> ret;
  for (auto &arr : ret) arr.resize(n+1);
  
  enum {t, z};
  auto &coefs(ret[z]); // GSL FFT works in-situ

  // all coefs zero...
  coefs=0.;
  // ... but the one defining the background ascent/descent
  coefs(1) = -.25 * z_hlf / si::metres; 
  // ... and the one defining the small-scale fluctuations
  if (freq != 0. * si::hertz) 
    coefs(int(2. * t_hlf * freq)) = -.5 * (ampl / si::metres);

  // assuring no sine components
  // (-cosine is used to assure symmetry in time)
  assert(
    max(coefs(blitz::Range(n/2-1, n-1))) == 0 && 
    min(coefs(blitz::Range(n/2-1, n-1))) == 0
  );

  // generating the signal with an inverse FFT
  int status = gsl_fft_halfcomplex_radix2_backward(coefs.data(), 1, n);
  assert(status == GSL_SUCCESS);

  for (int i=0; i < n+1; ++i)
    ret[t](i) = double(i) / n * 2 * t_hlf/si::seconds;

  ret[z] -= ret[z](0);
  ret[z](n) = 0;

  return ret;
/*
  // output
  {
    std::cout 
      << (double(i)/n * 2 * t_hlf/si::metres).value() 
      << " " 
      << arr(i)-arr(0) 
      << std::endl;
  }
  std::cout 
    << (2.*t_hlf/si::metres).value() 
    << " " 
    << 0 
    << std::endl;
*/
}
