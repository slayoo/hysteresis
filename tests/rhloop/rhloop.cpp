#include <string>
#include <stdexcept>
#include <set>

#include <blitz/array.h>
#include <gnuplot-iostream.h>

#include "../common/hdf5.hpp"

int main()
{
  Gnuplot gp;
  gp << "set encoding utf8\n";
  gp << "set term svg size 5.5*210,5.5*297 dynamic enhanced fsize 15\n";
  gp << "set output 'rhloop.svg'\n";
  gp << "set multiplot layout 3,3 columnsfirst\n";
  gp << "set grid xtics, ytics, mxtics, mytics\n";
  gp << "set mxtics 2\n";
  gp << "set mytics 2\n";

  for (auto &mlt : std::list<std::pair<double,double>>{{1., .1}, {1.,1.}, {1./2,1.}})
  {
    const auto 
      &mlt_r = mlt.first,
      &mlt_n = mlt.second;
    
    gp << "unset logscale\n";
    std::map<double, H5::H5File> cases;

    // simulations
    for (auto &mlt_t : std::set<double>{1.,2.,200}) 
    {
      std::ostringstream fname;
      fname << "refcase_t_x_" << mlt_t << "_rd_x_" << mlt_r << "_n_x_" << mlt_n << ".nc";
      {
	std::string cmd = "cp ../common/refcase.nc " + fname.str(); 
	system(cmd.c_str());
      }
      {
	H5::H5File file(fname.str().c_str(), H5F_ACC_RDWR);
	h5_setpar(file, "t_hlf", mlt_t * h5_getpar(file, "t_hlf"));
	h5_setpar(file, "rd",    mlt_r * h5_getpar(file, "rd"));
	h5_setpar(file, "N_stp", mlt_n * h5_getpar(file, "N_stp"));
      }
      {
	std::string cmd = "../../main " + fname.str(); 
	system(cmd.c_str());
      }
      cases[mlt_t] = H5::H5File(fname.str().c_str(), H5F_ACC_RDONLY);
    }

    double lw0 = 7., dlw = 1.05;
    {
      gp << "set key bottom right samplen .5\n";

      gp << "set xlabel 'S=RH-1 [%]' offset 0,.5\n";
      gp << "set ylabel 'displacement [m]' offset 1.5,0\n";
      gp << "set xtics .5\n";
      gp << "set ytics 50\n";
      gp << "set xrange [-.92:.92]\n";
      gp << "set yrange [0:160]\n";

      {
	std::ostringstream tmp;
	tmp << 1e6*h5_getpar(cases.begin()->second, "rd");
	gp << "set label 1 'r_d = " << tmp.str() << " µm' at graph .05,.85\n";
      }
      {
	std::ostringstream tmp;
	tmp << 1e-6*h5_getpar(cases.begin()->second, "N_stp");
	gp << "set label 2 'N_{STP} = " << tmp.str() << " cm^{-3}' at graph .05,.90\n";
      }
      gp << "set label 3 'ascent          ' at graph .05,.95\n";
      gp << "set label 4 '             / descent' at graph .05,.95 textcolor rgb 'orange'\n";

      gp << "plot 1./0 notitle";
      {
	int lw=lw0;
	for (auto &cs : cases)
	{
	  auto &h5 = cs.second;
	  double
	    t_hlf = h5_getpar(h5, "t_hlf"), 
	    z_hlf = h5_getpar(h5, "z_hlf");
	  std::ostringstream w_mean;
	  w_mean << (100 * z_hlf / t_hlf);
	  gp
	    << ",'-' using (100*($1-1)):($3<=" << 1.05*t_hlf << "?$2:1./0) with lines lc rgb 'black' lw " << lw << " title '<w>=" << w_mean.str() << " cm/s'"
	    << ",'-' using (100*($1-1)):($3>=" << t_hlf << "?$2:1./0) with lines lc rgb 'orange' lw " << lw/3. << " notitle";
	  lw-=dlw;
	}
      }
      gp << "\n";
      for (auto &cs : cases)
      {
	auto &h5 = cs.second;
	std::array<decltype(h5_getarr(h5, "t")), 3> data;
	data[0].reference(h5_getarr(h5, "RH"));
	data[1].reference(h5_getarr(h5, "z"));
	data[2].reference(h5_getarr(h5, "t"));
	for (int i=0; i<2; ++i) gp.send(data);
      }

      gp << "set ylabel 'S=RH-1 [%]' offset 1.5,0\n";
      gp << "set xlabel 'wet radius [µm]' offset 0,.5\n";
      gp << "set ytics .5\n";
      gp << "set xrange [.08:13.5]\n";
      gp << "set logscale x\n";
      gp << "set yrange [-.92:.92]\n";
      gp << "plot 1./0 notitle lc rgb 'black'";

      {
	int lw=lw0;
	for (auto &cs : cases)
	{ 
	  auto &h5 = cs.second;
	  double t_hlf = h5_getpar(h5, "t_hlf");
	  if (lw == lw0) 
	    gp
	      << ",'-' using ($1*1e6):(100*($4*$5-1)) with lines lc rgb 'grey' lw " << 1.75*lw0 << " title 'Köhler curve'"
	    ;
	  gp
	    << ",'-' using ($1*1e6):($3<=" << 1.05*t_hlf << "?100*($2-1):1./0) with lines lc rgb 'black' lw " << lw << " notitle"
	    << ",'-' using ($1*1e6):($3>=" << t_hlf << "?100*($2-1):1./0) with lines lc rgb 'orange'  lw " << lw/3. << " notitle"
	  ;
	  lw-=dlw;
	}
      }
      gp << "\n";

      bool wbg=true;
      for (auto &cs : cases)
      {
	auto &h5 = cs.second;
	std::array<decltype(h5_getarr(h5, "t")), 5> data;
	data[0].reference(h5_getarr(h5, "rw"));
	data[1].reference(h5_getarr(h5, "RH"));
	data[2].reference(h5_getarr(h5, "t"));
	data[3].reference(h5_getarr(h5, "kelvin"));
	data[4].reference(h5_getarr(h5, "raoult"));
	for (int i=0; i<(wbg?3:2); ++i) gp.send(data);
	wbg=false;
      }
    }

    gp << "set ylabel 'displacement [m]' offset 1.5,0\n";
    gp << "set yrange [25:75]\n";
    gp << "set ytics 25\n";
    gp << "plot 1./0 notitle lc rgb 'black'";
    {
      int lw=lw0;
      for (auto &cs : cases)
      { 
	auto &h5 = cs.second;
	double t_hlf = h5_getpar(h5, "t_hlf");
	gp
	  << ",'-' using ($1*1e6):($3<=" << 1.05*t_hlf << "?$2:1./0) with lines lc rgb 'black' lw " << lw << " notitle"
	  << ",'-' using ($1*1e6):($3>=" << t_hlf << "?$2:1./0) with lines lc rgb 'orange'  lw " << lw/3. << " notitle"
	;
	lw-=dlw;
      }
    }
    gp << "\n";

    bool wbg=true;
    for (auto &cs : cases)
    {
      auto &h5 = cs.second;
      std::array<decltype(h5_getarr(h5, "t")), 3> data;
      data[0].reference(h5_getarr(h5, "rw"));
      data[1].reference(h5_getarr(h5, "z"));
      data[2].reference(h5_getarr(h5, "t"));
      for (int i=0; i<2; ++i) gp.send(data);
    }
  }
}
