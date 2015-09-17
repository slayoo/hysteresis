#include <set>
#include <map>

#include <blitz/array.h>
#include <gnuplot-iostream.h> 

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include "../common/hdf5.hpp"

int main()
{
  std::array<std::string, 1> ixs{"p2"};
  const double min_RH = .99;

  // simulations 
  std::map<std::string,
    std::map<double, std::map<std::string, std::map<double, H5::H5File>>>
  > cases;

  char lll = 'a'-1;
  for (auto &mlt : std::list<std::pair<double,double>>{{1., .1}, {1.,1.}, {1./2,1.}})
  {                                                                             
    lll++;
    const auto                                                                  
      &mlt_r = mlt.first,                                                       
      &mlt_n = mlt.second;

    for (auto &mlt_t : std::list<double>{1.,2.,200})
    {
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
	  fname << "refcase_reltol_x_" << mlt_tol << "_ix=" << ix << "_t_x_" << mlt_t << "_r_x_" << mlt_r << "_n_x_" << mlt_n << ".nc";
	  {                                                                         
	    std::string cmd = "cp ../common/refcase.nc " + fname.str();             
	    system(cmd.c_str());                                                    
	  }                                                                         
	  double w;
          std::ostringstream label;
	  {                                                                         
	    H5::H5File file(fname.str().c_str(), H5F_ACC_RDWR);                     
	    h5_setpar(file, "reltol", mlt_tol * h5_getpar(file, "reltol"));
	    h5_setpar(file, "rd", mlt_r * h5_getpar(file, "rd"));
	    h5_setpar(file, "N_stp", mlt_n * h5_getpar(file, "N_stp"));
	    h5_setpar(file, "t_hlf", mlt_t * h5_getpar(file, "t_hlf"));
	    w = h5_getpar(file, "z_hlf") / h5_getpar(file, "t_hlf");
	    label 
              << "(" << lll << ") "
	      << "    "
	      << "N_{STP}=" << 1e-6 * h5_getpar(file, "N_stp") << " cm^{-3}"
	      << "    "
	      << "r_d=" << 1e6 * h5_getpar(file, "rd") << " µm";
	  }                                                                         
	  {                                                                         
	    std::string cmd = "../../main " + fname.str() + " " + ix;
	    system(cmd.c_str());                                                    
	  }                                                                         
	  cases[label.str()][w][ix][mlt_tol] = H5::H5File(fname.str().c_str(), H5F_ACC_RDONLY);   
	}
      }
    }
  }

  // analysis 
  std::map<std::string,
    std::map<double,
      std::map<std::string, std::array<blitz::Array<double, 1>,3>>
    >
  > plotdata;

  for (const auto &cs_lab : cases)
  {
    auto &lab = cs_lab.first;
    for (const auto &cs_w : cs_lab.second)
    {
      const auto &w = cs_w.first;
      for (const auto &ix : ixs)
      {
	for (auto &arr : plotdata[lab][w][ix]) 
	  arr.resize(cases[lab][w][ix].size());

	int csidx = 0;
	for (auto &cs : cases[lab][w][ix])
	{
	  auto &h5 = cs.second;
	  std::array<decltype(h5_getarr(h5, "t")), 2> data;
	  data[0].reference(h5_getarr(h5, "t"));
	  data[0] = blitz::where(h5_getarr(h5, "RH")>min_RH, h5_getarr(h5, "t"), NAN);
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

	  plotdata[lab][w][ix][0](csidx) = log10(h5_getpar(h5, "reltol"));
	  plotdata[lab][w][ix][1](csidx) = gsl_stats_quantile_from_sorted_data(dt_up.data(), 1, n_up, .1);
	  plotdata[lab][w][ix][2](csidx) = gsl_stats_quantile_from_sorted_data(dt_dn.data(), 1, n_dn, .1);
	  csidx++;
	}
      }
    }
  }

  // plotting
  Gnuplot gp;
  gp << "set term svg dynamic\n";
  gp << "set output 'dtrtol.svg'\n";
  gp << "set term svg size 5.5*210,5.5*297 dynamic enhanced fsize 15\n";        
  gp << "set multiplot layout 3,3\n";


  gp << "set label 2 'ascent          ' at graph .51,.05\n";                
  gp << "set label 3 '            / descent' at graph .51,.05 textcolor rgb 'orange'\n";

  gp << "set xlabel 'log_{10}(relative tolerance) [1]'\n";
  gp << "set ylabel '10^{th} percentile of timestep length distribution [s]' offset 3,1\n";
  gp << "set grid\n";
  gp << "set yrange [.0005:10]\n";
  gp << "set key top left reverse Left samplen .5\n";
  gp << "set logscale y\n";

  std::map<std::string, int> pts{{"id",1},{"ln",2},{"p2",4},{"p3",5}};
  std::map<std::string, double> lws{{"id",.5},{"ln",.5},{"p2",4},{"p3",.5}};
  std::map<std::string, std::string> tts{{"id","r"},{"ln","ln(r/1nm)"},{"p2","r^2"},{"p3","r^3"}};

  for (auto &cs_lab : plotdata)
  { 
    auto const &lab = cs_lab.first;
    gp << "set title '" << lab << "'\n";
    gp << "plot 1./0 not";
/*
    for (auto &ix : ixs)
      gp
	<< ",'-' using 1:2 with linespoints lc rgb " << (ix=="p2"?"'black'":"'gray'") << "  pt " << pts[ix] << " ps .666 lw " << lws[ix] << " title '" << tts[ix] << "'"
	<< ",'-' using 1:3 with linespoints lc rgb " << (ix=="p2"?"'orange'":"'gray'") << " pt " << pts[ix] << " ps .666 lw " << lws[ix] << " notitle '" << tts[ix] << "'"
     ;
*/
    double lw0=3., dlw=2.;
    double lw=lw0;
    for (auto &cs_w : plotdata[lab])
    {
      std::ostringstream w;
      w << 100*cs_w.first;
      gp 
	<< ",'-' using 1:2 with lines lc rgb 'black' lw " << lw << " title '<w>=" << w.str() << " cm/s'"
	<< ",'-' using 1:3 with lines lc rgb 'orange' lw " << lw << " notitle"
      ;
      lw+=dlw;
    }
    gp << "\n";
/*
    for (auto &ix : ixs)
      for (int i=0; i < 2; ++i) 
	gp.send(plotdata[lab][1.][ix]);
*/
    for (auto &cs_w : plotdata[lab])
      for (int i=0; i < 2; ++i) 
	gp.send(plotdata[lab][cs_w.first]["p2"]);
  }
}
