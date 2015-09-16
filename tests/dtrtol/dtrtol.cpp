#include <set>

#include <blitz/array.h>
#include <gnuplot-iostream.h> 

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include "../common/hdf5.hpp"

int main()
{
  std::array<std::string, 4> ixs{"id", "p2", "p3", "ln"};

  // simulations 
  std::map<std::string, std::map<double, H5::H5File>> cases;
  for (auto &ix : ixs)
  {  
    for (auto &mlt_tol : std::set<double>{
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
      fname << "refcase_reltol_x_" << mlt_tol << "_ix=" << ix << ".nc";           
      {                                                                         
	std::string cmd = "cp ../common/refcase.nc " + fname.str();             
	system(cmd.c_str());                                                    
      }                                                                         
      {                                                                         
	H5::H5File file(fname.str().c_str(), H5F_ACC_RDWR);                     
	h5_setpar(file, "reltol", mlt_tol * h5_getpar(file, "reltol"));
      }                                                                         
      {                                                                         
	std::string cmd = "../../main " + fname.str() + " " + ix;
	system(cmd.c_str());                                                    
      }                                                                         
      cases[ix][mlt_tol] = H5::H5File(fname.str().c_str(), H5F_ACC_RDONLY);   
    }
  }

  // analysis 
  std::map<std::string, std::array<blitz::Array<double, 1>,3>> plotdata;

  for (const auto &ix : ixs)
  {
    for (auto &arr : plotdata[ix]) 
      arr.resize(cases[ix].size());
    int csidx = 0;
    for (auto &cs : cases[ix])
    {
      auto &h5 = cs.second;
      std::array<decltype(h5_getarr(h5, "t")), 2> data;
      data[0].reference(h5_getarr(h5, "t"));
      data[0] = blitz::where(h5_getarr(h5, "RH")>.99, h5_getarr(h5, "t"), NAN);
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

      plotdata[ix][0](csidx) = log10(h5_getpar(h5, "reltol"));
      plotdata[ix][1](csidx) = gsl_stats_quantile_from_sorted_data(dt_up.data(), 1, n_up, .1);
      plotdata[ix][2](csidx) = gsl_stats_quantile_from_sorted_data(dt_dn.data(), 1, n_dn, .1);
      csidx++;
    }
  }

  // plotting
  Gnuplot gp;
  gp << "set term svg dynami\n";
  gp << "set output 'dtrtol.svg'\n";

  gp << "set label 2 'ascent          ' at graph .05,.92\n";                
  gp << "set label 3 '           / descent' at graph .05,.92 textcolor rgb 'orange'\n";

  gp << "set xlabel 'log_{10}(relative tolerance) [1]'\n";
  gp << "set ylabel '10^{th} percentile of timestep length distribution [s]'\n";
  gp << "set grid\n";
  gp << "set yrange [.005:.5]\n";
  gp << "set key bottom right samplen .5\n";
  gp << "set logscale y\n";

  std::map<std::string, int> pts{{"id",1},{"ln",2},{"p2",4},{"p3",5}};
  std::map<std::string, double> lws{{"id",.5},{"ln",3},{"p2",1},{"p3",1.5}};
  std::map<std::string, std::string> tts{{"id","r"},{"ln","ln(r/1µm)"},{"p2","r^2"},{"p3","r^3"}};

  gp << "plot 1./0 not";
  for (auto &ix : ixs)
    gp
      << ",'-' using 1:2 with linespoints lc rgb 'black'  pt " << pts[ix] << " ps .666 lw " << lws[ix] << " title '" << tts[ix] << "'"
      << ",'-' using 1:3 with linespoints lc rgb 'orange' pt " << pts[ix] << " ps .666 lw " << lws[ix] << " notitle '" << tts[ix] << "'"
   ;
  gp << "\n";
  for (auto &ix : ixs)
    for (int i=0; i < 2; ++i) 
      gp.send(plotdata[ix]);
}
