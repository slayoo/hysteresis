#include <set>

#include <blitz/array.h>
#include <gnuplot-iostream.h> 

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include "../common/hdf5.hpp"

int main()
{
  // simulations 
  std::map<double, H5::H5File> cases;
  for (auto &mlt_tol : std::set<double>{
    pow(5, 2),
    pow(5, 1),
    pow(5, 0),
    pow(5,-1),
    pow(5,-2),
    pow(5,-3),
    pow(5,-4),
    pow(5,-5),
    pow(5,-6),
    pow(5,-7),
    pow(5,-8),
  }) {
    std::ostringstream fname;                                                 
    fname << "refcase_reltol_x_" << mlt_tol << ".nc";           
    {                                                                         
      std::string cmd = "cp ../common/refcase.nc " + fname.str();             
      system(cmd.c_str());                                                    
    }                                                                         
    {                                                                         
      H5::H5File file(fname.str().c_str(), H5F_ACC_RDWR);                     
      h5_setpar(file, "reltol", mlt_tol * h5_getpar(file, "reltol"));
    }                                                                         
    {                                                                         
      std::string cmd = "../../main " + fname.str();                          
      system(cmd.c_str());                                                    
    }                                                                         
    cases[mlt_tol] = H5::H5File(fname.str().c_str(), H5F_ACC_RDONLY);   
  }

  // analysis 
  std::array<blitz::Array<double, 1>, 5> plotdata;
  for (auto &arr : plotdata) arr.resize(cases.size());
  int csidx = 0;
  for (auto &cs : cases)
  {
    auto &h5 = cs.second;
    std::array<decltype(h5_getarr(h5, "t")), 2> data;
    data[0].reference(h5_getarr(h5, "t"));
    data[0] = blitz::where(h5_getarr(h5, "RH")>.98, h5_getarr(h5, "t"), std::nanf(""));
    auto n = data[0].size();
    auto &dt = data[1];
    dt.resize(n-1);
    dt = data[0](blitz::Range(1,n-1)) - data[0](blitz::Range(0,n-2));
    data[0].reference(data[0](blitz::Range(1,n-1)));

    const double FLAG = 1e8;
    for (auto &val : dt) val = std::isfinite(val) ? val : FLAG; 

    blitz::Array<double,1> 
      dt_up(blitz::where(data[0] < h5_getpar(h5, "t_hlf"), dt, FLAG)),
      dt_dn(blitz::where(data[0] > h5_getpar(h5, "t_hlf"), dt, FLAG));

    gsl_sort(dt_up.data(), 1, dt_up.size());
    gsl_sort(dt_dn.data(), 1, dt_dn.size());

    for (auto &val : dt_up) val = val == FLAG ? NAN : val; 
    for (auto &val : dt_dn) val = val == FLAG ? NAN : val; 

    auto 
      n_up = gsl_stats_min_index(dt_up.data(), 1, n),
      n_dn = gsl_stats_min_index(dt_dn.data(), 1, n);

    plotdata[0](csidx) = log10(h5_getpar(h5, "reltol"));
    plotdata[1](csidx) = gsl_stats_quantile_from_sorted_data(dt_up.data(), 1, n_up, .05);
    plotdata[2](csidx) = gsl_stats_quantile_from_sorted_data(dt_up.data(), 1, n_up, .5);
    plotdata[3](csidx) = gsl_stats_quantile_from_sorted_data(dt_dn.data(), 1, n_dn, .05);
    plotdata[4](csidx) = gsl_stats_quantile_from_sorted_data(dt_dn.data(), 1, n_dn, .5);
    csidx++;
  }

  // plotting
  Gnuplot gp;
  gp << "set term svg dynami\n";
  gp << "set output 'dtrtol.svg'\n";

  gp << "set xlabel 'log_{10}(relative tolerance) [1]'\n";
  gp << "set ylabel 'timestep [s]'\n";
  gp << "set grid\n";
  gp << "set key bottom right samplen .5\n";
  gp << "set logscale y\n";

  gp << "plot "
     << "  '-' using 1:2 with linesp lc rgb 'black' lw 3 title 'ascent / 5^{th} percentile',"
     << "  '-' using 1:3 with linesp lc rgb 'black'  title 'ascent / median',"
     << "  '-' using 1:4 with linesp lc rgb 'orange' lw 3 title 'descent / 5^{th} percentile',"
     << "  '-' using 1:5 with linesp lc rgb 'orange' title 'descent / median'"
     << "\n";
  for (int i=0; i < 4; ++i) gp.send(plotdata);
}
