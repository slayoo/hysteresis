#include <array>
#include <gsl/gsl_fft_halfcomplex.h>
#include <blitz/array.h>
#include <boost/units/systems/si.hpp>

#include "gslerr.hpp"

struct invfft_t
{
  struct params_t
  {
    const boost::units::quantity<boost::units::si::length>    z_hlf;
    const boost::units::quantity<boost::units::si::time>      t_hlf;
    const boost::units::quantity<boost::units::si::frequency> freq;
    const boost::units::quantity<boost::units::si::length>    ampl;
  };
  
  std::array<blitz::Array<double, 1>, 2> operator()(const params_t &params) const
  {
    gslerr_init();

    namespace si = boost::units::si;

    // at least 64 points per cycle
    const int n = pow(2, ceil(log2(64. * (
      params.freq == 0. * si::hertz 
	? 1.
	: (params.freq * (2. * params.t_hlf)).value()
    ))));

    std::array<blitz::Array<double, 1>, 2> ret;
    for (auto &arr : ret) arr.resize(n+1);
    
    enum {t, z};
    auto &coefs(ret[z]); // GSL FFT works in-situ

    // all coefs zero...
    coefs=0.;
    // ... but the one defining the background ascent/descent
    coefs(1) = -.25 * params.z_hlf / si::metres; 
    // ... and the one defining the small-scale fluctuations
    if (params.freq != 0. * si::hertz) 
      coefs(int(2. * params.t_hlf * params.freq)) = -.5 * (params.ampl / si::metres);

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
      ret[t](i) = double(i) / n * 2 * params.t_hlf/si::seconds;

    ret[z] -= ret[z](0);
    ret[z](n) = 0;

    return ret;
  }
};
