#include <stdexcept>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <blitz/array.h>
#include <boost/units/systems/si.hpp>

using namespace boost::units;

const auto z_hlf = 100. * si::metres;
const auto t_hlf = 100. * si::seconds;
const auto freq = .46 * si::hertz;
const auto ampl = 1. * si::metres;

void gslerr2xcpt(const char* reason, const char* file, int line, int gsl_errno)
{
  throw std::runtime_error(std::string(file) + ": " + std::string(reason));
}

int main()
{
  // at least 64 points per cycle
  const int n = pow(2, ceil(log2(64. * freq * (2. * t_hlf))));

  blitz::Array<double, 1> arr(n);

  // all coefs zero...
  arr=0;
  // ... but the one defining the background ascent/descent
  arr(1) = -.25 * z_hlf / si::metres; 
  // ... and the one defining the small-scale fluctuations
  arr(int(2. * t_hlf * freq)) = -.5 * (ampl / si::metres);

  // assuring no sine components
  // (-cosine is used to assure symmetry in time)
  assert(
    max(arr(blitz::Range(n/2-1, n))) == 0 && 
    min(arr(blitz::Range(n/2-1, n))) == 0
  );

  // generating the signal with an inverse FFT
  gsl_set_error_handler(gslerr2xcpt);
  int status = gsl_fft_halfcomplex_radix2_backward(
    arr.data(), 
    arr.stride(0), 
    arr.extent(0)
  );
  assert(status == GSL_SUCCESS);

  // output
  for (int i=0; i < n; ++i)
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
}
