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
#pragma once
#include "../core/common/types.hpp"
#include <fstream>
#include <iostream>
#include <netcdf>

#ifdef USE_MPI
#include <netcdf_par.h>
#include <mpi.h>
#endif

using namespace netCDF;
using namespace netCDF::exceptions;

namespace io_muphys {

#ifdef __SINGLE_PRECISION
#define NC_REAL_TYPE NC_FLOAT
#define NC_PUT_VARA(ncid,varid,start,count,ptr) nc_put_vara_float(ncid,varid,start,count,ptr)
#define NC_GET_VARA(ncid,varid,start,count,ptr) nc_get_vara_float(ncid,varid,start,count,ptr)
using NCreal_t = NcFloat;
#else
#define NC_REAL_TYPE NC_DOUBLE
#define NC_PUT_VARA(ncid,varid,start,count,ptr) nc_put_vara_double(ncid,varid,start,count,ptr)
#define NC_GET_VARA(ncid,varid,start,count,ptr) nc_get_vara_double(ncid,varid,start,count,ptr)
using NCreal_t = NcDouble;
#endif

void parse_args(std::string &file, string &outfile, size_t &itime, real_t &dt, real_t &qnc,
                int argc, char **argv);

void input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                  const std::string input, size_t &ncells, size_t &nlev,
                  size_t itime);
void input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                  const std::string input, size_t &ncells, size_t &nlev);

void output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                   const std::string output, array_1d_t<real_t> &v,
                   size_t &ncells, size_t &nlev, int &deflate_level);
void output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                   const std::string output, std::map<std::string, NcVarAtt>,
                   array_1d_t<real_t> &v, size_t &ncells, size_t &nlev,
                   int &deflate_level);

void read_fields(const std::string input_file, size_t &itime, size_t &ncells,
                 size_t &nlev, array_1d_t<real_t> &z, array_1d_t<real_t> &t,
                 array_1d_t<real_t> &p, array_1d_t<real_t> &rho,
                 array_1d_t<real_t> &qv, array_1d_t<real_t> &qc,
                 array_1d_t<real_t> &qi, array_1d_t<real_t> &qr,
                 array_1d_t<real_t> &qs, array_1d_t<real_t> &qg);

void write_fields(const string output_file, size_t &ncells, size_t &nlev,
                  array_1d_t<real_t> &t, array_1d_t<real_t> &qv,
                  array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                  array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                  array_1d_t<real_t> &qg, array_1d_t<real_t> &prr_gsp,
                  array_1d_t<real_t> &pri_gsp, array_1d_t<real_t> &prs_gsp,
                  array_1d_t<real_t> &prg_gsp, array_1d_t<real_t> &pre_gsp,
                  array_1d_t<real_t> &pflx);
void write_fields(const string output_file, const string input_file,
                  size_t &ncells, size_t &nlev, array_1d_t<real_t> &t,
                  array_1d_t<real_t> &qv, array_1d_t<real_t> &qc,
                  array_1d_t<real_t> &qi, array_1d_t<real_t> &qr,
                  array_1d_t<real_t> &qs, array_1d_t<real_t> &qg,
                  array_1d_t<real_t> &prr_gsp, array_1d_t<real_t> &pri_gsp,
                  array_1d_t<real_t> &prs_gsp, array_1d_t<real_t> &prg_gsp,
                  array_1d_t<real_t> &pre_gsp, array_1d_t<real_t> &pflx);
} // namespace io_muphys

#ifdef USE_MPI
namespace io_muphys {
  #ifdef __SINGLE_PRECISION
  #define NC_REAL_TYPE NC_FLOAT
  #define NC_PUT_VARA(ncid,varid,start,count,ptr) nc_put_vara_float(ncid,varid,start,count,ptr)
  #define NC_GET_VARA(ncid,varid,start,count,ptr) nc_get_vara_float(ncid,varid,start,count,ptr)
  #else
  #define NC_REAL_TYPE NC_DOUBLE
  #define NC_PUT_VARA(ncid,varid,start,count,ptr) nc_put_vara_double(ncid,varid,start,count,ptr)
  #define NC_GET_VARA(ncid,varid,start,count,ptr) nc_get_vara_double(ncid,varid,start,count,ptr)
  #endif
  void parse_args_mpi_rank0(string &file, string &outfile, size_t &itime, real_t &dt, real_t &qnc, int argc, char **argv);
  void parse_args_mpi(string &file, string &outfile, size_t &itime, real_t &dt, real_t &qnc, int argc, char **argv);
  // Read a vector variable (level and cell dimensions) at given time index.
  void input_vector_mpi(const std::string &filename, array_1d_t<real_t> &v,
                        const std::string &var_name, size_t ncells,
                        size_t nlev, size_t itime, MPI_Comm comm, MPI_Info info);
  
  // Read all fields (wrapper for multiple input_vector_mpi calls).
  void read_fields_mpi(const string input_file, size_t &itime,
                        size_t &ncells, size_t &nlev, array_1d_t<real_t> &z,
                        array_1d_t<real_t> &t, array_1d_t<real_t> &p,
                        array_1d_t<real_t> &rho, array_1d_t<real_t> &qv,
                        array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                        array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                        array_1d_t<real_t> &qg, MPI_Comm comm = MPI_COMM_WORLD, MPI_Info info = MPI_INFO_NULL);

  void output_vector_par(int ncid, int varid, size_t itime, size_t start_cell, size_t ncell_loc, size_t nlev, const array_1d_t<real_t>& arr);

  void output_vector_par(int ncid, int varid, size_t start_cell, size_t ncell_loc, size_t nlev, const array_1d_t<real_t>& arr);

  void write_fields_mpi(const std::string &output_file, size_t ncells, size_t nlev,
                        const array_1d_t<real_t> &t, const array_1d_t<real_t> &qv,
                        const array_1d_t<real_t> &qc, const array_1d_t<real_t> &qi,
                        const array_1d_t<real_t> &qr, const array_1d_t<real_t> &qs,
                        const array_1d_t<real_t> &qg, const array_1d_t<real_t> &prr_gsp,
                        const array_1d_t<real_t> &pri_gsp, const array_1d_t<real_t> &prs_gsp,
                        const array_1d_t<real_t> &prg_gsp, const array_1d_t<real_t> &pre_gsp,
                        const array_1d_t<real_t> &pflx, int deflate_level = 0,
                        MPI_Comm comm = MPI_COMM_WORLD, MPI_Info info = MPI_INFO_NULL);
} // namespace io_muphys
#endif