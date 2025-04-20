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

void capture_by_col(array_1d_t<real_t>& var, array_1d_t<real_t>& buffer, size_t rows, size_t cols, size_t j) {
   for (size_t i = 0; i < rows; ++i) {
      size_t oned_vec_index = i * cols + j;
      buffer[i] = var[oned_vec_index];
   }
}

void capture_by_col_single(array_1d_t<real_t>& var, array_1d_t<real_t>& buffer, size_t j) {
   buffer[0] = var[j];
}

void write_by_col_single(array_1d_t<real_t>& var, array_1d_t<real_t>& buffer, size_t j) {
   var[j] = buffer[0];
}

void write_by_col(array_1d_t<real_t>& var, array_1d_t<real_t>& buffer, size_t rows, size_t cols, size_t j) {
   for (size_t i = 0; i < rows; ++i) {
      size_t oned_vec_index = i * cols + j;
      var[oned_vec_index] = buffer[i];
   }
}

int main(int argc, char *argv[]) {
   // MPI_Proc
   MPI_Init(&argc, &argv);
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   const int root = 0;

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

   size_t ncols = ivend - ivbeg;
   array_2d_t<real_t> buffer(11, array_1d_t<real_t>(nlev * nvec, ZERO));
   array_2d_t<real_t> gsp(5, array_1d_t<real_t>(ncells, ZERO));
   for (size_t ii = 0; ii < multirun; ++ii){
      for (size_t jj = ivbeg; jj < ivend; ++jj) {
         capture_by_col(dz, buffer[0], nlev, nvec, jj);
         capture_by_col(t, buffer[1], nlev, nvec, jj);
         capture_by_col(rho, buffer[2], nlev, nvec, jj);
         capture_by_col(p, buffer[3], nlev, nvec, jj);
         capture_by_col(qv, buffer[4], nlev, nvec, jj);
         capture_by_col(qc, buffer[5], nlev, nvec, jj);
         capture_by_col(qi, buffer[6], nlev, nvec, jj);
         capture_by_col(qr, buffer[7], nlev, nvec, jj);
         capture_by_col(qs, buffer[8], nlev, nvec, jj);
         capture_by_col(qg, buffer[9], nlev, nvec, jj);
         capture_by_col(pflx, buffer[10], nlev, nvec, jj);

         capture_by_col_single(prr_gsp, gsp[0], jj);
         capture_by_col_single(pri_gsp, gsp[1], jj);
         capture_by_col_single(prs_gsp, gsp[2], jj);
         capture_by_col_single(prg_gsp, gsp[3], jj);
         capture_by_col_single(pre_gsp, gsp[4], jj);

         size_t temp_nvec  = 1;
         size_t temp_ivbeg = 0;
         size_t temp_ivend = 1;

         graupel(temp_nvec, kend, temp_ivbeg, temp_ivend, kbeg, dt, buffer[0], buffer[1], buffer[2], buffer[3],
                 buffer[4], buffer[5], buffer[6], buffer[7], buffer[8], buffer[9], qnc_1,
                  gsp[0], gsp[1],gsp[2], gsp[3], gsp[4], buffer[10]);

         write_by_col_single(prr_gsp, gsp[0], jj);
         write_by_col_single(pri_gsp, gsp[1], jj);
         write_by_col_single(prs_gsp, gsp[2], jj);
         write_by_col_single(prg_gsp, gsp[3], jj);
         write_by_col_single(pre_gsp, gsp[4], jj);

         write_by_col(dz, buffer[0], nlev, nvec, jj);
         write_by_col(t, buffer[1], nlev, nvec, jj);
         write_by_col(rho, buffer[2], nlev, nvec, jj);
         write_by_col(p, buffer[3], nlev, nvec, jj);
         write_by_col(qv, buffer[4], nlev, nvec, jj);
         write_by_col(qc, buffer[5], nlev, nvec, jj);
         write_by_col(qi, buffer[6], nlev, nvec, jj);
         write_by_col(qr, buffer[7], nlev, nvec, jj);
         write_by_col(qs, buffer[8], nlev, nvec, jj);
         write_by_col(qg, buffer[9], nlev, nvec, jj);
         write_by_col(pflx, buffer[10], nlev, nvec, jj);
       }
   } 
   auto end_time = std::chrono::steady_clock::now();
   auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time);

   io_muphys::write_fields(output_file, ncells, nlev, t, qv, qc, qi, qr, qs,
                              qg, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp, pflx);

   std::cout << "time taken : " << duration.count() << " milliseconds"
            << std::endl;

   MPI_Finalize();
   return 0;
}
