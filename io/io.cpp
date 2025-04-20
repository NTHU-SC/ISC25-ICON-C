// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------
//
#include "io.hpp"
#include <algorithm>
#include <map>

#define NC_CHECK(status)                                                    \
  if ((status) != NC_NOERR) {                                              \
    std::cerr << "NetCDF Error: " << nc_strerror(status) << std::endl;    \
    MPI_Abort(MPI_COMM_WORLD, status);                                     \
}

static const int NC_ERR = 2;
static const std::string BASE_VAR = "zg";

void io_muphys::parse_args(string &file, string &outfile, size_t &itime, real_t &dt, real_t &qnc,
                           int argc, char **argv) {
  file = "aes-new-gr_moderate-dt30s_atm_3d_ml_20080801T000000Z.nc";
  itime = 0; /* default to the first timestep in the input file */
  dt = 30.0;
  qnc = 100.0;
  char *end = nullptr;

  if (argc > 1) {
    file = argv[1];
  }
  //cout << "input file: " << file << "\n";

  if (argc > 2) {
    outfile = argv[2];
  }
  //cout << "output file: " << outfile << "\n";

  if (argc > 3) {
    itime = stoi(argv[2]);
  }
  //cout << "itime: " << itime << "\n";

  if (argc > 4) {
    dt = type_converter<real_t>(argv[3], &end);
  }
  //cout << "dt: " << dt << "\n";

  if (argc > 5) {
    qnc = type_converter<real_t>(argv[4], &end);
  }
  //cout << "qnc: " << qnc << endl;
}

/* read-in time-constant data fields without a time dimension */
void io_muphys::input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                             const string input, size_t &ncells, size_t &nlev) {
  NcVar var;
  v.resize(ncells * nlev);
  /*  access the input variable */
  try {
    var = datafile.getVar(input);
  } catch (NcNotVar &e) {
    cout << "FAILURE in accessing " << input << " (no time dimension) *******"
         << endl;
    e.what();
    e.errorCode();
  }
  /*  read-in input field values */
  try {
    array_1d_t<size_t> startp = {0, 0};
    array_1d_t<size_t> count = {nlev, ncells};
    var.getVar(startp, count, v.data());
  } catch (NcNotVar &e) {
    cout << "FAILURE in reading values from " << input
         << " (no time dimensions) *******" << endl;
    e.what();
    e.errorCode();
  }
}

void io_muphys::input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                             const string input, size_t &ncells, size_t &nlev,
                             size_t itime) {
  NcVar att = datafile.getVar(input);
  try {
    v.resize(ncells * nlev);
    if (att.isNull()) {
      throw NC_ERR;
    }
    array_1d_t<size_t> startp = {itime, 0, 0};
    array_1d_t<size_t> count = {1, nlev, ncells};
    att.getVar(startp, count, v.data());
  } catch (NcException &e) {
    e.what();
    cout << "FAILURE in reading " << input << " **************************"
         << endl;
    throw NC_ERR;
  }
}

void io_muphys::output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                              const string output, array_1d_t<real_t> &v,
                              size_t &ncells, size_t &nlev,
                              int &deflate_level) {
  // fortran:column major while c++ is row major
  NCreal_t ncreal_t;
  netCDF::NcVar var = datafile.addVar(output, ncreal_t, dims);

  for (size_t i = 0; i < nlev; ++i) {
    var.putVar({i, 0}, {1, ncells}, &v[i * ncells]);
  }
  if (deflate_level > 0) {
    var.setCompression(true, false, deflate_level);
  }
}

void io_muphys::output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                              const string output,
                              std::map<std::string, NcVarAtt> varAttributes,
                              array_1d_t<real_t> &v, size_t &ncells,
                              size_t &nlev, int &deflate_level) {
  // fortran:column major while c++ is row major
  NCreal_t ncreal_t;
  netCDF::NcVar var = datafile.addVar(output, ncreal_t, dims);

  for (size_t i = 0; i < nlev; ++i) {
    var.putVar({i, 0}, {1, ncells}, &v[i * ncells]);
  }
  if (deflate_level > 0) {
    var.setCompression(true, false, deflate_level);
  }
  /* Add given attribues to the output variables (string, only) */
  for (auto &attribute_name : {"standard_name", "long_name", "units",
                               "coordinates", "CDI_grid_type"}) {
    auto attribute = varAttributes[attribute_name];

    std::string dataValues = "default";
    attribute.getValues(dataValues);

    var.putAtt(attribute.getName(), attribute.getType(),
               attribute.getAttLength(), dataValues.c_str());
  }
}

void io_muphys::read_fields(const string input_file, size_t &itime,
                            size_t &ncells, size_t &nlev, array_1d_t<real_t> &z,
                            array_1d_t<real_t> &t, array_1d_t<real_t> &p,
                            array_1d_t<real_t> &rho, array_1d_t<real_t> &qv,
                            array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                            array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                            array_1d_t<real_t> &qg) {

  NcFile datafile(input_file, NcFile::read);

  /*  read in the dimensions from the base variable: zg
   *  1st) vertical
   *  2nd) horizontal
   */
  auto baseDims = datafile.getVar(BASE_VAR).getDims();
  nlev = baseDims[0].getSize();
  ncells = baseDims[1].getSize();

  io_muphys::input_vector(datafile, z, "zg", ncells, nlev);
  io_muphys::input_vector(datafile, t, "ta", ncells, nlev, itime);

  io_muphys::input_vector(datafile, p, "pfull", ncells, nlev, itime);
  io_muphys::input_vector(datafile, rho, "rho", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qv, "hus", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qc, "clw", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qi, "cli", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qr, "qr", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qs, "qs", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qg, "qg", ncells, nlev, itime);

  datafile.close();
}
void io_muphys::write_fields(
    string output_file, size_t &ncells, size_t &nlev, array_1d_t<real_t> &t,
    array_1d_t<real_t> &qv, array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
    array_1d_t<real_t> &qr, array_1d_t<real_t> &qs, array_1d_t<real_t> &qg,
    array_1d_t<real_t> &prr_gsp, array_1d_t<real_t> &pri_gsp,
    array_1d_t<real_t> &prs_gsp, array_1d_t<real_t> &prg_gsp,
    array_1d_t<real_t> &pre_gsp, array_1d_t<real_t> &pflx) {
  NcFile datafile(output_file, NcFile::replace);
  NcDim ncells_dim = datafile.addDim("ncells", ncells);
  NcDim nlev_dim = datafile.addDim("height", nlev);
  std::vector<NcDim> dims = {nlev_dim, ncells_dim};
  int deflate_level = 0;
  size_t onelev = 1;
  NcDim onelev_dim = datafile.addDim("height1", onelev);
  std::vector<NcDim> dims1d = {onelev_dim, ncells_dim};

  io_muphys::output_vector(datafile, dims, "ta", t, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "hus", qv, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "clw", qc, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "cli", qi, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qr", qr, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qs", qs, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qg", qg, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "pflx", pflx, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "prr_gsp", prr_gsp, ncells, onelev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "prs_gsp", prs_gsp, ncells, onelev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "pri_gsp", pri_gsp, ncells, onelev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "prg_gsp", prg_gsp, ncells, onelev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "pre_gsp", pre_gsp, ncells, onelev,
                           deflate_level);

  datafile.close();
}

[[maybe_unused]] static void copy_coordinate_variables_if_present(NcFile &datafile, NcFile &inputfile,
                                     std::vector<std::string> coordinates) {
  for (auto &coordinate_name : coordinates) {
    auto coordinate = inputfile.getVar(coordinate_name);

    /* skip if variable wasn't found*/
    if (coordinate.isNull())
      continue;

    /* copy possible new dimensions from input coordinates to the output
     * befor adding the related data variables */
    for (NcDim &dim : coordinate.getDims()) {
      auto currentDims = datafile.getDims();
      /* map.contains() would be better, but requires c++20 */
      if (auto search = currentDims.find(dim.getName());
          search == currentDims.end())
        datafile.addDim(dim.getName(), dim.getSize());
    }
    auto var = datafile.addVar(coordinate.getName(), coordinate.getType(),
                               coordinate.getDims());

    for (auto &[attribute_name, attribute] : coordinate.getAtts()) {
      std::string dataValues = "default";
      attribute.getValues(dataValues);
      var.putAtt(attribute.getName(), attribute.getType(),
                 attribute.getAttLength(), dataValues.c_str());
    }

    /* Calculate the total size of the varibale */
    auto dimensions = coordinate.getDims();
    struct Prod {
      void operator()(NcDim n) { prod *= n.getSize(); }
      size_t prod{1};
    };
    size_t totalSize =
        std::for_each(dimensions.cbegin(), dimensions.cend(), Prod()).prod;

    /* Create a one-dimensional vector to store its values */
    std::vector<real_t> oneDimensionalVariable(totalSize);

    /* Copy the original values to the output */
    coordinate.getVar(oneDimensionalVariable.data());
    var.putVar(oneDimensionalVariable.data());
  }
}

void io_muphys::write_fields(
    string output_file, string input_file, size_t &ncells, size_t &nlev,
    array_1d_t<real_t> &t, array_1d_t<real_t> &qv, array_1d_t<real_t> &qc,
    array_1d_t<real_t> &qi, array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
    array_1d_t<real_t> &qg, array_1d_t<real_t> &prr_gsp,
    array_1d_t<real_t> &pri_gsp, array_1d_t<real_t> &prs_gsp,
    array_1d_t<real_t> &prg_gsp, array_1d_t<real_t> &pre_gsp,
    array_1d_t<real_t> &pflx) {
  NcFile datafile(output_file, NcFile::replace);
  NcFile inputfile(input_file, NcFile::read);
  auto baseDims = inputfile.getVar(BASE_VAR).getDims();
  NcDim nlev_dim =
      datafile.addDim(baseDims[0].getName(), baseDims[0].getSize());
  NcDim ncells_dim =
      datafile.addDim(baseDims[1].getName(), baseDims[1].getSize());

  copy_coordinate_variables_if_present(datafile, inputfile,
                                       {"clon", "clon_bnds", "clat",
                                        "clat_bnds", baseDims[0].getName(),
                                        "height_bnds"});
  /*TODO  height_bnds might have a different name */

  std::vector<NcDim> dims = {nlev_dim, ncells_dim};
  int deflate_level = 0;
  size_t onelev = 1;
  NcDim onelev_dim = datafile.addDim("height1", onelev);
  std::vector<NcDim> dims1d = {onelev_dim, ncells_dim};

  io_muphys::output_vector(datafile, dims, "ta",
                           inputfile.getVar("ta").getAtts(), t, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "hus",
                           inputfile.getVar("hus").getAtts(), qv, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "clw",
                           inputfile.getVar("clw").getAtts(), qc, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "cli",
                           inputfile.getVar("cli").getAtts(), qi, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qr",
                           inputfile.getVar("qr").getAtts(), qr, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qs",
                           inputfile.getVar("qs").getAtts(), qs, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qg",
                           inputfile.getVar("qg").getAtts(), qg, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "pflx", pflx, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "prr_gsp", prr_gsp, ncells, onelev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "prs_gsp", prs_gsp, ncells, onelev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "pri_gsp", pri_gsp, ncells, onelev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "prg_gsp", prg_gsp, ncells, onelev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims1d, "pre_gsp", pre_gsp, ncells, onelev,
                           deflate_level);
  inputfile.close();
  datafile.close();
}

[[maybe_unused]] static void copy_coordinate_variables(NcFile &datafile, NcFile &inputfile,
                                      array_1d_t<std::string> coordinates) {
  for (auto &coordinate_name : coordinates) {
    auto coordinate = inputfile.getVar(coordinate_name);
    /* copy possible new dimensions from input coordinates to the output
     * befor adding the related data variables */
    for (NcDim &dim : coordinate.getDims()) {
      auto currentDims = datafile.getDims();
      /* map.contains() would be better, but requires c++20 */
      if (auto search = currentDims.find(dim.getName());
          search == currentDims.end())
        datafile.addDim(dim.getName(), dim.getSize());
    }
    auto var = datafile.addVar(coordinate.getName(), coordinate.getType(),
                               coordinate.getDims());

    for (auto &[attribute_name, attribute] : coordinate.getAtts()) {
      std::string dataValues = "default";
      attribute.getValues(dataValues);
      var.putAtt(attribute.getName(), attribute.getType(),
                 attribute.getAttLength(), dataValues.c_str());
    }

    /* Calculate the total size of the varibale */
    auto dimensions = coordinate.getDims();
    struct Prod {
      void operator()(NcDim n) { prod *= n.getSize(); }
      size_t prod{1};
    };
    size_t totalSize =
        std::for_each(dimensions.cbegin(), dimensions.cend(), Prod()).prod;

    /* Create a one-dimensional vector to store its values */
    array_1d_t<real_t> oneDimensionalVariable(totalSize);

    /* Copy the original values to the output */
    coordinate.getVar(oneDimensionalVariable.data());
    var.putVar(oneDimensionalVariable.data());
  }
}

namespace io_muphys {

  void input_scalar_mpi(const std::string &filename, real_t &value,
                        const std::string &var_name, size_t itime,
                        MPI_Comm comm) {
    int ncid, varid;
    NC_CHECK(nc_open_par(filename.c_str(), NC_NOWRITE | NC_MPIIO,
                         comm, MPI_INFO_NULL, &ncid));
    NC_CHECK(nc_inq_varid(ncid, var_name.c_str(), &varid));
    NC_CHECK(nc_var_par_access(ncid, varid, NC_COLLECTIVE));
    size_t start[1] = {itime};
    size_t count[1] = {1};
    NC_CHECK(nc_get_vara_double(ncid, varid, start, count, &value));
    NC_CHECK(nc_close(ncid));
  }
  
  void input_vector_mpi(const std::string &filename, array_1d_t<real_t> &v,
                        const std::string &var_name, size_t ncells,
                        size_t nlev, size_t itime, MPI_Comm comm) {
    int ncid, varid;
    int mpi_rank, mpi_size;
    MPI_Comm_rank(comm, &mpi_rank);
    MPI_Comm_size(comm, &mpi_size);
    NC_CHECK(nc_open_par(filename.c_str(), NC_NOWRITE | NC_MPIIO,
                         comm, MPI_INFO_NULL, &ncid));
    NC_CHECK(nc_inq_varid(ncid, var_name.c_str(), &varid));
    NC_CHECK(nc_var_par_access(ncid, varid, NC_COLLECTIVE));
  
    size_t base = ncells / mpi_size;
    size_t rem = ncells % mpi_size;
    size_t local_ncells = base + (mpi_rank < rem ? 1 : 0);
    size_t start_cell = mpi_rank * base + std::min<size_t>(mpi_rank, rem);
  
    size_t start[3] = {itime, 0, start_cell};
    size_t count[3] = {1, nlev, local_ncells};
    v.resize(nlev * local_ncells);
    NC_CHECK(nc_get_vara_double(ncid, varid, start, count, v.data()));
    NC_CHECK(nc_close(ncid));
  }
  
  void read_fields_mpi(const std::string &input_file, size_t itime,
                       size_t ncells, size_t nlev,
                       array_1d_t<real_t> &z, array_1d_t<real_t> &t,
                       array_1d_t<real_t> &p, array_1d_t<real_t> &rho,
                       array_1d_t<real_t> &qv, array_1d_t<real_t> &qc,
                       array_1d_t<real_t> &qi, array_1d_t<real_t> &qr,
                       array_1d_t<real_t> &qs, array_1d_t<real_t> &qg,
                       array_1d_t<real_t> &prr_gsp,
                       array_1d_t<real_t> &pri_gsp,
                       array_1d_t<real_t> &prs_gsp,
                       array_1d_t<real_t> &prg_gsp,
                       array_1d_t<real_t> &pre_gsp,
                       array_1d_t<real_t> &pflx,
                       MPI_Comm comm) {
    input_vector_mpi(input_file, z, "zg", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, t, "ta", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, p, "pfull", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, rho, "rho", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, qv, "hus", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, qc, "clw", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, qi, "cli", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, qr, "qr", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, qs, "qs", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, qg, "qg", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, prr_gsp, "prr_gsp", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, pri_gsp, "pri_gsp", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, prs_gsp, "prs_gsp", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, prg_gsp, "prg_gsp", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, pre_gsp, "pre_gsp", ncells, nlev, itime, comm);
    input_vector_mpi(input_file, pflx, "pflx", ncells, nlev, itime, comm);
  }
  
  void write_fields_mpi(const std::string &output_file,
                        size_t ncells, size_t nlev,
                        const array_1d_t<real_t> &t,
                        const array_1d_t<real_t> &qv,
                        const array_1d_t<real_t> &qc,
                        const array_1d_t<real_t> &qi,
                        const array_1d_t<real_t> &qr,
                        const array_1d_t<real_t> &qs,
                        const array_1d_t<real_t> &qg,
                        const array_1d_t<real_t> &prr_gsp,
                        const array_1d_t<real_t> &pri_gsp,
                        const array_1d_t<real_t> &prs_gsp,
                        const array_1d_t<real_t> &prg_gsp,
                        const array_1d_t<real_t> &pre_gsp,
                        const array_1d_t<real_t> &pflx,
                        MPI_Comm comm) {
    int ncid;
    int dim_time, dim_level, dim_cell;
    int mpi_rank, mpi_size;
    MPI_Comm_rank(comm, &mpi_rank);
    MPI_Comm_size(comm, &mpi_size);
  
    NC_CHECK(nc_create_par(output_file.c_str(), NC_NETCDF4 | NC_MPIIO | NC_CLOBBER,
                           comm, MPI_INFO_NULL, &ncid));
  
    NC_CHECK(nc_def_dim(ncid, "time", NC_UNLIMITED, &dim_time));
    NC_CHECK(nc_def_dim(ncid, "level", nlev, &dim_level));
    NC_CHECK(nc_def_dim(ncid, "cell", ncells, &dim_cell));
  
    int dims[3] = {dim_time, dim_level, dim_cell};
    std::vector<std::pair<const char*, const array_1d_t<real_t>&>> vars = {
      {"ta", t}, {"hus", qv}, {"clw", qc}, {"cli", qi}, {"qr", qr}, {"qs", qs}, {"qg", qg},
      {"prr_gsp", prr_gsp}, {"pri_gsp", pri_gsp}, {"prs_gsp", prs_gsp},
      {"prg_gsp", prg_gsp}, {"pre_gsp", pre_gsp}, {"pflx", pflx}
    };
  
    std::vector<int> var_ids(vars.size());
  
    for (size_t i = 0; i < vars.size(); ++i) {
      NC_CHECK(nc_def_var(ncid, vars[i].first, NC_DOUBLE, 3, dims, &var_ids[i]));
      NC_CHECK(nc_var_par_access(ncid, var_ids[i], NC_COLLECTIVE));
    }
  
    NC_CHECK(nc_enddef(ncid));
  
    size_t base = ncells / mpi_size;
    size_t rem = ncells % mpi_size;
    size_t local_ncells = base + (mpi_rank < rem ? 1 : 0);
    size_t start_cell = mpi_rank * base + std::min<size_t>(mpi_rank, rem);
    size_t start[3] = {0, 0, start_cell};
    size_t count[3] = {1, nlev, local_ncells};
  
    for (size_t i = 0; i < vars.size(); ++i) {
      const auto& var = vars[i].second;
      NC_CHECK(nc_put_vara_double(ncid, var_ids[i], start, count, var.data()));
    }
  
    NC_CHECK(nc_close(ncid));
  }
  
} // namespace io_muphys