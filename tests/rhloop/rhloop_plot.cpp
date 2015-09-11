#include <string>
#include <stdexcept>

#include <H5Cpp.h>
#include <blitz/array.h>
#include <gnuplot-iostream.h>

auto h5load(const H5::H5File &file, const std::string &var)
{
  auto dset = file.openDataSet(var);
  auto dspc = dset.getSpace();

  const int rank = 1;
  assert(dspc.getSimpleExtentNdims() == rank);

  hsize_t n;
  dspc.getSimpleExtentDims(&n, NULL);
  blitz::Array<float, 1> data(n);
  dset.read(data.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(rank, &n), dspc);

  return data;
}

int main()
{
  H5::H5File file("rhloop.nc", H5F_ACC_RDONLY);
  Gnuplot gp;

  gp << "set term svg dynamic enhanced mouse standalone fsize 18\n";
  gp << "set output 'rhloop.svg'\n";
  gp << "set xlabel 'RH [Pa/Pa]'\n";
  gp << "set ylabel 'z [m]'\n";
  gp << "set xtics rotate by -65\n";
  gp << "set grid\n";
  gp << "set arrow from 1, graph 0 to 1, graph 1 nohead\n";

  std::array<decltype(h5load(file, "RH")),2> data;
  data[0].reference(h5load(file, "RH"));
  data[1].reference(h5load(file, "z"));

  gp << "plot '-' with lines notitle\n";
  gp.send(data);
}
