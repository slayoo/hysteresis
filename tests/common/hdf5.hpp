#pragma once

#include <string>

#include <H5Cpp.h>
#include <blitz/array.h>

double h5_getpar(
  const H5::H5File &file, 
  const std::string &name
) {
  double value;
  auto attr = file.openGroup("/").openAttribute(name);
  attr.read(attr.getDataType(), &value);
  return value;
}

void h5_setpar(
  const H5::H5File &file, 
  const std::string &name, 
  const double &value
) {
  auto attr = file.openGroup("/").openAttribute(name);
  attr.write(attr.getDataType(), &value);
}

auto h5_getarr(const H5::H5File &file, const std::string &var)                     
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
