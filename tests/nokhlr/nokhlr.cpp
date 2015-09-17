#include <set>

#include <blitz/array.h>
#include <gnuplot-iostream.h>

#include "../common/hdf5.hpp"

int main()
{
  // simulations                                                                
  std::map<double, std::map<bool, H5::H5File>> cases;     
  for (auto &mlt_t : std::set<double>{.5,1,250})
  {
    for (auto &khl : std::set<bool>{0,1}) 
    { 
      std::ostringstream fname;                                                 
      fname << "refcase_khl=" << khl << "_t_x_" << mlt_t << ".nc";         
      {                                                                         
	std::string cmd = "cp ../common/refcase.nc " + fname.str();             
	system(cmd.c_str());                                                    
      }                                                                         
      {                                                                         
	H5::H5File file(fname.str().c_str(), H5F_ACC_RDWR);                     
	h5_setpar(file, "kelvin", double(khl));
	h5_setpar(file, "raoult", double(khl));
	h5_setpar(file, "RH0", 1.-1e-6);
	h5_setpar(file, "reltol", 1e-8);
	h5_setpar(file, "n_cycl", 1.);
	h5_setpar(file, "t_hlf", mlt_t * h5_getpar(file, "t_hlf"));
      }                                                                         
      {                                                                         
	std::string cmd = "../../main " + fname.str();
	system(cmd.c_str());                                                    
      }                                                                         
      cases[mlt_t][khl] = H5::H5File(fname.str().c_str(), H5F_ACC_RDONLY);     
    }
  }

  // plotting
  Gnuplot gp;
  gp << "set term svg dynamic\n";
  gp << "set output 'nokhlr.svg'\n";
  gp << "set grid\n";
  gp << "set key above samplen .5\n";
  gp << "set ylabel 'displacement [z]'\n";
  gp << "set yrange [-10:160]\n";
  gp << "set xlabel 'S=RH-1 [%]'\n";


  gp << "plot 1./0 not";
  int lw = 3;
  for (auto &c_w : cases)
  {
    for (auto &cs : c_w.second)
    {
      auto &h5 = cs.second;
      bool khl = cs.first;

      std::ostringstream label;
      label << "Köhler: " << (khl?"on ":"off");

      double                                                                
	t_hlf = h5_getpar(h5, "t_hlf"),                                     
	z_hlf = h5_getpar(h5, "z_hlf");                                     
      std::ostringstream w_mean;                                            
      w_mean << (100 * z_hlf / t_hlf); 
      gp << ", '-' using (100*($1-1)):2 with lines lw " << lw << " lc rgb " << (khl?"'red'":"'blue'") << " title '" << label.str() << " (w=" << w_mean.str() << " cm/s)'";
    }
    lw -= 1;
  }
  gp << "\n";

  for (auto &c_w : cases)
  {
    for (auto &cs : c_w.second)
    {
      auto &h5 = cs.second;
      std::array<blitz::Array<float,1>, 2> data;
      data[0].reference(h5_getarr(h5, "RH"));
      data[1].reference(h5_getarr(h5, "z"));
      gp.send(data);
    }
  }
}
