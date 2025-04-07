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
#include <chrono>
#include <cstdlib>
#include <iostream>

#include "core/common/graupel.hpp"
#include "core/common/types.hpp"
#include "core/common/utils.hpp"
#include "io/io.hpp"
#include <chrono>

int main(int argc, char *argv[]) {
  // Parameters from the command line
  string file;
  string output_file;
  size_t itime;
  real_t dt, qnc, qnc_1;
  io_muphys::parse_args(file, output_file, itime, dt, qnc, argc, argv);

  // Parameters from the input file
  size_t ncells, nlev;
  array_1d_t<real_t> z, t, p, rho, qv, qc, qi, qr, qs, qg;

  // Pre-calculated parameters
  array_1d_t<real_t> dz;
  // Extra fields required to call graupel
  array_1d_t<real_t> pflx, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp;
  // start-end indices
  size_t kend, kbeg, ivend, ivbeg, nvec;

  const string input_file = file;
  io_muphys::read_fields(input_file, itime, ncells, nlev, z, t, p, rho, qv, qc,
                         qi, qr, qs, qg);
  utils_muphys::calc_dz(z, dz, ncells, nlev);

  prr_gsp.resize(ncells, ZERO);
  pri_gsp.resize(ncells, ZERO);
  prs_gsp.resize(ncells, ZERO);
  prg_gsp.resize(ncells, ZERO);
  pre_gsp.resize(ncells, ZERO);
  pflx.resize(ncells * nlev, ZERO);

  kbeg = 0;
  kend = nlev;
  ivbeg = 0;
  ivend = ncells;
  nvec = ncells;
  qnc_1 = qnc;

  size_t multirun = 0;

  if (std::getenv("MULTI_GRAUPEL")){
     multirun = atoi(std::getenv("MULTI_GRAUPEL"));
     }
  else {
     multirun = 1;
     }
  std::cout << "multirun =" << multirun << std::endl;

  auto start_time = std::chrono::steady_clock::now();

  for (size_t ii = 0; ii < multirun; ++ii){
    graupel(nvec, kend, ivbeg, ivend, kbeg, dt, dz, t, rho, p, qv, qc, qi, qr, qs,
            qg, qnc_1, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp, pflx);
 } 
  auto end_time = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time);

  io_muphys::write_fields(output_file, ncells, nlev, t, qv, qc, qi, qr, qs,
                            qg, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp, pflx);

  std::cout << "time taken : " << duration.count() << " milliseconds"
            << std::endl;


  return 0;
}
