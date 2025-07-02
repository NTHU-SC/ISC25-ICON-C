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
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[]) {
   //auto start_time = std::chrono::steady_clock::now();
   // Mpi parameters initilization
   MPI_Init(&argc, &argv);
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   // Parameters from the command line
   string file;
   string output_file;
   size_t itime;
   real_t dt, qnc, qnc_1;
   
   if (!rank)
      io_muphys::parse_args_mpi_rank0(file, output_file, itime, dt, qnc, argc, argv);
   else
      io_muphys::parse_args_mpi(file, output_file, itime, dt, qnc, argc, argv);

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

   io_muphys::read_fields_mpi(input_file, itime, ncells, nlev, z, t, p, rho, qv, qc, qi, qr, qs, qg);

   size_t base = ncells / size;
   size_t rem  = ncells % size;
   size_t ncell_loc = base + (rank < rem ? 1 : 0);
   size_t start_cell = rank * base + std::min(rank, (int)rem);

   utils_muphys::calc_dz(z, dz, ncell_loc, nlev);

   prr_gsp.resize(ncell_loc, ZERO);
   pri_gsp.resize(ncell_loc, ZERO);
   prs_gsp.resize(ncell_loc, ZERO);
   prg_gsp.resize(ncell_loc, ZERO);
   pre_gsp.resize(ncell_loc, ZERO);
   pflx.resize(ncell_loc * nlev, ZERO);

   kbeg = 0;
   kend = nlev;
   ivbeg = 0;
   ivend = ncell_loc;
   nvec = ncell_loc;
   qnc_1 = qnc;
   

   size_t multirun = 0;

   if (std::getenv("MULTI_GRAUPEL")){
      multirun = atoi(std::getenv("MULTI_GRAUPEL"));
   }
   else {
      multirun = 1;
   }

   if (!rank)
      std::cout << "multirun =" << multirun << std::endl;

   auto start_time = std::chrono::steady_clock::now();
   for (size_t ii = 0; ii < multirun; ++ii){
      graupel(nvec, kend, ivbeg, ivend, kbeg, dt, dz, t, rho, p, qv, qc, qi, qr, qs,
              qg, qnc_1, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp, pflx);
   }                     
   auto end_time = std::chrono::steady_clock::now();
   io_muphys::write_fields_mpi(output_file, ncells, nlev, t, qv, qc, qi, qr, qs,
         qg, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp, pflx);
   
   //auto end_time = std::chrono::steady_clock::now();
   auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
             end_time - start_time);
   if (!rank)
      std::cout << "time taken : " << duration.count() << " milliseconds" << std::endl;

   MPI_Finalize();
   return 0;
}