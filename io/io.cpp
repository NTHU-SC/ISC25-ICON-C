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
  cout << "input file: " << file << "\n";

  if (argc > 2) {
    outfile = argv[2];
  }
  cout << "output file: " << outfile << "\n";

  if (argc > 3) {
    itime = stoi(argv[2]);
  }
  cout << "itime: " << itime << "\n";

  if (argc > 4) {
    dt = type_converter<real_t>(argv[3], &end);
  }
  cout << "dt: " << dt << "\n";

  if (argc > 5) {
    qnc = type_converter<real_t>(argv[4], &end);
  }
  cout << "qnc: " << qnc << endl;
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
#ifdef USE_MPI
namespace io_muphys {
  void parse_args_mpi_rank0(string &file, string &outfile, size_t &itime, real_t &dt, real_t &qnc,
                           int argc, char **argv) {
    file = "aes-new-gr_moderate-dt30s_atm_3d_ml_20080801T000000Z.nc";
    itime = 0; /* default to the first timestep in the input file */
    dt = 30.0;
    qnc = 100.0;
    char *end = nullptr;

    if (argc > 1) {
      file = argv[1];
    }
    cout << "input file: " << file << "\n";

    if (argc > 2) {
      outfile = argv[2];
    }
    cout << "output file: " << outfile << "\n";

    if (argc > 3) {
      itime = stoi(argv[2]);
    }
    cout << "itime: " << itime << "\n";

    if (argc > 4) {
      dt = type_converter<real_t>(argv[3], &end);
    }
    cout << "dt: " << dt << "\n";

    if (argc > 5) {
      qnc = type_converter<real_t>(argv[4], &end);
    }
    cout << "qnc: " << qnc << endl;
  }

  void parse_args_mpi(string &file, string &outfile, size_t &itime, real_t &dt, real_t &qnc,
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
  // Helper: read a 3D variable [time, lev, cell] in parallel
  void input_vector_mpi(int ncid,
                        const char *name,
                        size_t itime,
                        size_t start_cell,
                        size_t ncell_loc,
                        size_t nlev,
                        array_1d_t<real_t> &arr) {
      int varid;
      // find variable ID
      if (nc_inq_varid(ncid, name, &varid)) 
          throw std::runtime_error(std::string("Variable not found: ") + name);
      // set collective I/O
      nc_var_par_access(ncid, varid, NC_INDEPENDENT);
      // define the hyperslab: [time, level, cell]
      size_t start[3] = { itime, 0, start_cell };
      size_t count[3] = { 1, nlev, ncell_loc };
      // allocate local buffer
      arr.resize(nlev * ncell_loc);
      // read
      if (NC_GET_VARA(ncid, varid, start, count, arr.data()))
          throw std::runtime_error(std::string("Failed to read var: ") + name);
  }

  // Helper: read a 2D variable [lev, cell] (static) in parallel
  void input_vector_mpi(int ncid,
                        const char *name,
                        size_t start_cell,
                        size_t ncell_loc,
                        size_t nlev,
                        array_1d_t<real_t> &arr) {
      int varid;
      // find variable ID
      if (nc_inq_varid(ncid, name, &varid)) 
          throw std::runtime_error(std::string("Variable not found: ") + name);
      // independent I/O is fine for static fields
      nc_var_par_access(ncid, varid, NC_INDEPENDENT);
      // define hyperslab: [level, cell]
      size_t start[2] = { 0, start_cell };
      size_t count[2] = { nlev, ncell_loc };
      // allocate local buffer
      arr.resize(nlev * ncell_loc);
      // read
      if (NC_GET_VARA(ncid, varid, start, count, arr.data()))
          throw std::runtime_error(std::string("Failed to read var: ") + name);
  }
  // Helper: write a 3D variable [time, level, cell] in parallel
  void output_vector_par(int ncid,
                         int varid,
                         size_t itime,
                         size_t start_cell,
                         size_t ncell_loc,
                         size_t nlev,
                         const array_1d_t<real_t>& arr) {

      // define hyperslab: [time, level, cell]
      size_t startp[3] = { itime, 0, start_cell };
      size_t countp[3] = { 1, nlev, ncell_loc };

      // write data (assuming arr is ordered as [level][cell])
      if (NC_PUT_VARA(ncid, varid, startp, countp, arr.data())) {
          throw std::runtime_error("Failed to write var: " + std::to_string(varid));
      }
  }
  // Helper: write a 2D variable [level, cell] in parallel
  void output_vector_par(int ncid,
                         int varid,
                         size_t start_cell,
                         size_t ncell_loc,
                         size_t nlev,
                         const array_1d_t<real_t>& arr) {
      // write data block (assuming arr is ordered as [level][cell])
      for (size_t lvl = 0; lvl < nlev; ++lvl) {
          size_t startp[2] = { lvl, start_cell };
          size_t countp[2] = { 1, ncell_loc };
          const real_t* data_ptr = arr.data() + lvl * ncell_loc;
          if (NC_PUT_VARA(ncid, varid, startp, countp, data_ptr)) {
              throw std::runtime_error("Failed to write var: " + std::to_string(varid));
          }
      }
  }
  void read_fields_mpi(const string input_file, size_t &itime,
                            size_t &ncells, size_t &nlev, array_1d_t<real_t> &z,
                            array_1d_t<real_t> &t, array_1d_t<real_t> &p,
                            array_1d_t<real_t> &rho, array_1d_t<real_t> &qv,
                            array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                            array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                            array_1d_t<real_t> &qg, MPI_Comm comm, MPI_Info info) {
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int ncid;
    // open file in parallel mode
    if (nc_open_par(input_file.c_str(), NC_NOWRITE | NC_MPIIO, comm, info, &ncid)) {
        throw std::runtime_error("Failed to open NetCDF file in parallel");
    }

    // inquire dimensions from "zg": [lev, cell]
    
    int varid_zg;
    nc_inq_varid(ncid, BASE_VAR.c_str(), &varid_zg);
    size_t dimlen_lev, dimlen_cell;
    int dimids[2];
    nc_inq_vardimid(ncid, varid_zg, dimids);  // Gets dimension IDs for zg
    nc_inq_dimlen(ncid, dimids[0], &dimlen_lev); 
    nc_inq_dimlen(ncid, dimids[1], &dimlen_cell);

    nlev   = dimlen_lev;
    ncells = dimlen_cell;

    // compute local cell block
    size_t base = ncells / nprocs;
    size_t rem  = ncells % nprocs;
    size_t ncell_loc = base + (rank < rem ? 1 : 0);
    size_t start_cell = rank * base + std::min(rank, (int)rem);

    // read static vertical grid
    io_muphys::input_vector_mpi(ncid, "zg", start_cell, ncell_loc, nlev, z);

    // read time‐dependent fields
    io_muphys::input_vector_mpi(ncid, "ta",  itime, start_cell, ncell_loc, nlev, t);
    io_muphys::input_vector_mpi(ncid, "pfull",itime, start_cell, ncell_loc, nlev, p);
    io_muphys::input_vector_mpi(ncid, "rho", itime, start_cell, ncell_loc, nlev, rho);
    io_muphys::input_vector_mpi(ncid, "hus", itime, start_cell, ncell_loc, nlev, qv);
    io_muphys::input_vector_mpi(ncid, "clw", itime, start_cell, ncell_loc, nlev, qc);
    io_muphys::input_vector_mpi(ncid, "cli", itime, start_cell, ncell_loc, nlev, qi);
    io_muphys::input_vector_mpi(ncid, "qr",  itime, start_cell, ncell_loc, nlev, qr);
    io_muphys::input_vector_mpi(ncid, "qs",  itime, start_cell, ncell_loc, nlev, qs);
    io_muphys::input_vector_mpi(ncid, "qg",  itime, start_cell, ncell_loc, nlev, qg);

    // close file
    nc_close(ncid);
  }

  // Write 3D and 2D fields in parallel, splitting horizontal dimension across ranks
  void write_fields_mpi(const std::string &output_file,
                        size_t ncells,
                        size_t nlev,
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
                        int deflate_level,
                        MPI_Comm comm,
                        MPI_Info info) {
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    // Compute local block for cell dimension
    size_t base = ncells / nprocs;
    size_t rem  = ncells % nprocs;
    size_t ncell_loc = base + (rank < (int)rem ? 1 : 0);
    size_t start_cell = rank * base + std::min<size_t>(rank, rem);

    int ncid;
    // create file in parallel mode
    if (nc_create_par(output_file.c_str(),
                      NC_NETCDF4 | NC_CLOBBER | NC_MPIIO | NC_CLASSIC_MODEL,
                      comm, info, &ncid)) {
        throw std::runtime_error("Failed to create NetCDF file in parallel");
    }

    // define dimensions: height and cell
    int dimid_height, dimid_cell;

    if (nc_def_dim(ncid, "ncells", ncells, &dimid_cell)) {
      throw std::runtime_error("Failed to define cell dimension");
  }

    if (nc_def_dim(ncid, "height", nlev, &dimid_height)) {
        throw std::runtime_error("Failed to define height dimension");
    }

    // define 1D height dimension for surface fields
    int dimid_height1;
    size_t onelev = 1;
    if (nc_def_dim(ncid, "height1", onelev, &dimid_height1)) {
        throw std::runtime_error("Failed to define height1 dimension");
    }

    std::map<std::string, int> varids;

    auto def_var = [&](const char* name, int dim0, int dim1) {
      int varid;
      int dimids[2] = {dim0, dim1};
      if (nc_def_var(ncid, name, NC_REAL_TYPE, 2, dimids, &varid))
          throw std::runtime_error("define failed");
      if (deflate_level > 0)
          nc_def_var_deflate(ncid, varid, 0, 1, deflate_level);
      nc_var_par_access(ncid, varid, NC_INDEPENDENT);
      varids[name] = varid;
    };

    def_var("ta", dimid_height, dimid_cell);
    def_var("hus", dimid_height, dimid_cell);
    def_var("clw", dimid_height, dimid_cell);
    def_var("cli", dimid_height, dimid_cell);
    def_var("qr", dimid_height, dimid_cell);
    def_var("qs", dimid_height, dimid_cell);
    def_var("qg", dimid_height, dimid_cell);
    def_var("pflx", dimid_height, dimid_cell);

    def_var("prr_gsp", dimid_height1, dimid_cell);
    def_var("prs_gsp", dimid_height1, dimid_cell);
    def_var("pri_gsp", dimid_height1, dimid_cell);
    def_var("prg_gsp", dimid_height1, dimid_cell);
    def_var("pre_gsp", dimid_height1, dimid_cell);

    // exit define mode
    nc_enddef(ncid);

    // 2D fields [height x cell]
    io_muphys::output_vector_par(ncid, varids["ta"], start_cell, ncell_loc, nlev, t);
    io_muphys::output_vector_par(ncid, varids["hus"], start_cell, ncell_loc, nlev, qv);
    io_muphys::output_vector_par(ncid, varids["clw"], start_cell, ncell_loc, nlev, qc);
    io_muphys::output_vector_par(ncid, varids["cli"], start_cell, ncell_loc, nlev, qi);
    io_muphys::output_vector_par(ncid, varids["qr"], start_cell, ncell_loc, nlev, qr);
    io_muphys::output_vector_par(ncid, varids["qs"], start_cell, ncell_loc, nlev, qs);
    io_muphys::output_vector_par(ncid, varids["qg"], start_cell, ncell_loc, nlev, qg);
    io_muphys::output_vector_par(ncid, varids["pflx"], start_cell, ncell_loc, nlev, pflx);

    // surface fields [height1 x cell]
    io_muphys::output_vector_par(ncid, varids["prr_gsp"], start_cell, ncell_loc, onelev, prr_gsp);
    io_muphys::output_vector_par(ncid, varids["prs_gsp"], start_cell, ncell_loc, onelev, prs_gsp);
    io_muphys::output_vector_par(ncid, varids["pri_gsp"], start_cell, ncell_loc, onelev, pri_gsp);
    io_muphys::output_vector_par(ncid, varids["prg_gsp"], start_cell, ncell_loc, onelev, prg_gsp);
    io_muphys::output_vector_par(ncid, varids["pre_gsp"], start_cell, ncell_loc, onelev, pre_gsp);

    // close file
    nc_close(ncid);
  }
} // namespace io_muphys
#endif