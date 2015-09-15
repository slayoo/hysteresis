#pragma once

#include <string>
#include <set>

#include <H5Cpp.h>

// the C HDF5 API (for dimension scales)
#include "hdf5_hl.h"

class hdf5io_t
{
  public:

  enum { 
    ix_t, ix_z, ix_RH, ix_rw, ix_T, ix_kelvin, ix_raoult,
    n_vars = 7
  };

  private:

  H5::H5File file;
  std::array<H5::DataSet, n_vars> dsets;
  const int n_dims = 1; // only time
  hsize_t shape = 1;

  public:

  auto getpar(const std::string &key) const
  {
    auto attr = file.openGroup("/").openAttribute(key);
    double data;
    attr.read(attr.getDataType(), &data);
    return data;
  }

  void putrec(const std::array<double, n_vars> &record)
  {

    const hsize_t 
      count = 1,
      offst = shape - 1;
    const H5::DataSpace dspc_data(n_dims, &count);

    for (int ix=0; ix < n_vars; ++ix)
    {
      dsets[ix].extend(&shape);
      auto dspc_trgt = dsets[ix].getSpace();
      dspc_trgt.selectHyperslab(H5S_SELECT_SET, &count, &offst);
      dsets[ix].write(&record[ix], H5::PredType::NATIVE_DOUBLE, dspc_data, dspc_trgt);
    }
    shape += 1;
  }

  // ctor
  hdf5io_t(const std::string &filename)
    : file(filename, H5F_ACC_RDWR)
  {
    // TODO? H5::Exception::dontPrint(); 

    // defining the file structure
    hsize_t 
      limit = H5S_UNLIMITED,
      chunk = 1;
    
    H5::DSetCreatPropList params;
    params.setChunk(1, &chunk); // required for unlimited dimension
    params.setDeflate(1);

    for (const auto &ixnm : std::map<int, const char*>({
      {ix_z , "z" },
      {ix_t , "t" },
      {ix_RH, "RH"},
      {ix_rw, "rw"},
      {ix_T,  "T"},
      {ix_kelvin, "kelvin"},
      {ix_raoult, "raoult"}
    })) {
      dsets[ixnm.first] = file.createDataSet(
        ixnm.second, 
        H5::PredType::NATIVE_FLOAT, // that's float, not double - just to save space
        H5::DataSpace(1, &shape, &limit), 
        params
      );
    }

    // phony_dim_0 -> time
    herr_t status;
    status = H5DSset_scale(dsets[ix_t].getId(), "time");
    assert(status == 0);
  }
};
