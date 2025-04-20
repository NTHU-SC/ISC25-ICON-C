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

   if (!rank)
      std::cout << "multirun =" << multirun << std::endl;

   // Step 1: Compute number of columns assigned to each rank
   array_1d_t<int> col_counts(size), col_displs(size);
   int base = ncells / size, rem = ncells % size;

   for (int r = 0; r < size; ++r) 
      col_counts[r] = base + (r < rem? 1: 0);

   col_displs[0] = 0;
   for (int r = 1; r < size; ++r)
      col_displs[r] = col_displs[r - 1] + col_counts[r - 1];
   
   int local_cols = col_counts[rank];

   MPI_Request reqs[16 * size];
   array_2d_t<real_t> buffer(11, array_1d_t<real_t>(nlev * local_cols));
   array_2d_t<real_t> gsp(5, array_1d_t<real_t>(local_cols));

   int req_count = 0;

   auto start_time = std::chrono::steady_clock::now();

   // Step 2: Scatter data
   if (size != 1 && !rank) {
      int sizes[2] = {nlev, ncells};
      for (int r = 1; r < size; ++r) {
         int start_col = col_displs[r];
         int cols = col_counts[r];
         
         int subsizes[2] = {nlev, cols};
         int starts[2] = {0, start_col};

         MPI_Datatype subarray;
         MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_REAL_T, &subarray);
         MPI_Type_commit(&subarray);

         MPI_Isend(dz.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(t.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(rho.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(p.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(qv.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(qc.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(qi.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(qr.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(qs.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(qg.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(pflx.data(), 1, subarray, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);

         MPI_Isend(prr_gsp.data() + start_col, cols, MPI_REAL_T, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(pri_gsp.data() + start_col, cols, MPI_REAL_T, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(prs_gsp.data() + start_col, cols, MPI_REAL_T, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(prg_gsp.data() + start_col, cols, MPI_REAL_T, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Isend(pre_gsp.data() + start_col, cols, MPI_REAL_T, r, 0, MPI_COMM_WORLD, &reqs[req_count++]);

         MPI_Type_free(&subarray);
      }
      
      for (size_t i = 0; i < nlev; ++i) {
         for (size_t j = 0; j < local_cols; ++j) {
            buffer[0][i * local_cols + j] = dz[i * nvec + j];
            buffer[1][i * local_cols + j] = t[i * nvec + j];
            buffer[2][i * local_cols + j] = rho[i * nvec + j];
            buffer[3][i * local_cols + j] = p[i * nvec + j];
            buffer[4][i * local_cols + j] = qv[i * nvec + j];
            buffer[5][i * local_cols + j] = qc[i * nvec + j];
            buffer[6][i * local_cols + j] = qi[i * nvec + j];
            buffer[7][i * local_cols + j] = qr[i * nvec + j];
            buffer[8][i * local_cols + j] = qs[i * nvec + j];
            buffer[9][i * local_cols + j] = qg[i * nvec + j];
            buffer[10][i * local_cols + j] = pflx[i * nvec + j];
         }
      }
      for (size_t i = 0; i < local_cols; ++i) {
         gsp[0][i] = prr_gsp[i];
         gsp[1][i] = pri_gsp[i];
         gsp[2][i] = prs_gsp[i];
         gsp[3][i] = prg_gsp[i];
         gsp[4][i] = pre_gsp[i];
      }
      MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
   } else if (size != 1){
      for (int i = 0; i < 11; ++i)
         MPI_Irecv(buffer[i].data(), nlev * local_cols, MPI_REAL_T, 0, 0, MPI_COMM_WORLD, &reqs[req_count++]);
      
      for (int i = 0; i < 5; ++i)
         MPI_Irecv(gsp[i].data(), local_cols, MPI_REAL_T, 0, 0, MPI_COMM_WORLD, &reqs[req_count++]);

      MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
   }
 
   size_t local_nvec = col_counts[rank];
   size_t local_ivbeg = 0;
   size_t local_ivend = local_nvec;

   // Step 3: Calculate result in each rank

   if (size == 1) {
      for (size_t ii = 0; ii < multirun; ++ii){
         graupel(nvec, kend, ivbeg, ivend, kbeg, dt, dz, t, rho, p, qv, qc, qi, qr, qs,
                 qg, qnc_1, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp, pflx);
      } 
   } else {
      for (size_t ii = 0; ii < multirun; ++ii){
         graupel(local_nvec, kend, local_ivbeg, local_ivend, kbeg, dt, buffer[0], buffer[1], buffer[2],
                 buffer[3], buffer[4], buffer[5], buffer[6], buffer[7], buffer[8],
                 buffer[9], qnc_1, gsp[0], gsp[1], gsp[2], gsp[3], gsp[4], buffer[10]);
      } 
   }

   req_count = 0;

   // Step 4: Gather
   if (size != 1 && !rank) {
      // gather data
      int sizes[2] = {nlev, ncells};
      for (int r = 1; r < size; ++r) {
         int start_col = col_displs[r];
         int cols = col_counts[r];

         int subsizes[2] = {nlev, cols};
         int starts[2] = {0, start_col};

         MPI_Datatype subarray;
         MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_REAL_T, &subarray);
         MPI_Type_commit(&subarray);

         MPI_Irecv(dz.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(t.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(rho.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(p.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(qv.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(qc.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(qi.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(qr.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(qs.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(qg.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(pflx.data(), 1, subarray, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);

         MPI_Irecv(prr_gsp.data() + start_col, cols, MPI_REAL_T, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(pri_gsp.data() + start_col, cols, MPI_REAL_T, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(prs_gsp.data() + start_col, cols, MPI_REAL_T, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(prg_gsp.data() + start_col, cols, MPI_REAL_T, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);
         MPI_Irecv(pre_gsp.data() + start_col, cols, MPI_REAL_T, r, 1, MPI_COMM_WORLD, &reqs[req_count++]);

         MPI_Type_free(&subarray);
      }

      for (size_t i = 0; i < nlev; ++i) {
         for (size_t j = 0; j < local_cols; ++j) {
            dz[i * nvec + j]     = buffer[0][i * local_cols + j];
            t[i * nvec + j]      = buffer[1][i * local_cols + j];
            rho[i * nvec + j]    = buffer[2][i * local_cols + j];
            p[i * nvec + j]      = buffer[3][i * local_cols + j];
            qv[i * nvec + j]     = buffer[4][i * local_cols + j];
            qc[i * nvec + j]     = buffer[5][i * local_cols + j];
            qi[i * nvec + j]     = buffer[6][i * local_cols + j];
            qr[i * nvec + j]     = buffer[7][i * local_cols + j];
            qs[i * nvec + j]     = buffer[8][i * local_cols + j];
            qg[i * nvec + j]     = buffer[9][i * local_cols + j];
            pflx[i * nvec + j]   = buffer[10][i * local_cols + j];
         }
      }
      for (size_t i = 0; i < local_cols; ++i) {
         prr_gsp[i] = gsp[0][i];
         pri_gsp[i] = gsp[1][i];
         prs_gsp[i] = gsp[2][i];
         prg_gsp[i] = gsp[3][i];
         pre_gsp[i] = gsp[4][i];
      }
      
      MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
   } else if (size != 1) {
      for (int i = 0; i < 11; ++i)
         MPI_Isend(buffer[i].data(), nlev * local_cols, MPI_REAL_T, 0, 1, MPI_COMM_WORLD, &reqs[req_count++]);
   
      for (int i = 0; i < 5; ++i)
         MPI_Isend(gsp[i].data(), local_cols, MPI_REAL_T, 0, 1, MPI_COMM_WORLD, &reqs[req_count++]);

      MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
   } 

   // Step 5: Output
   if (!rank) {
      // write fields
      io_muphys::write_fields(output_file, ncells, nlev, t, qv, qc, qi, qr, qs,
         qg, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pre_gsp, pflx);
   }

   auto end_time = std::chrono::steady_clock::now();
   auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
       end_time - start_time);

   if (!rank)
      std::cout << "time taken : " << duration.count() << " milliseconds" << std::endl;

   MPI_Finalize();
   return 0;
}
