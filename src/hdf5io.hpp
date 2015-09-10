#pragma once

#include <string>

#include <H5Cpp.h>

class hdf5io
{
  std::unique_ptr<H5::H5File> hdfp;

  public:

  

  // ctor
  hdf5io(const std::string &file)
  {
    hdfp.reset(new H5::H5File(file, H5F_ACC_RDWR));
  }
};
