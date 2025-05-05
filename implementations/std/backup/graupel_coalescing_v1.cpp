
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
#include "core/common/graupel.hpp"
#include <iostream>
#include <numeric>
#include <execution>
#include <algorithm>
#include <array>

using namespace property;
using namespace thermo;
using namespace transition;
using namespace idx;
using namespace graupel_ct;

#define T_SIZE 1024
#define ceil(a, b) ((a + b - 1) / b)

inline constexpr std::array<size_t, 6> make_qx_ind_array() {
  return {
      idx::qx_ind[0],
      idx::qx_ind[1],
      idx::qx_ind[2],
      idx::qx_ind[3],
      idx::qx_ind[4],
      idx::qx_ind[5]
  };
}

inline constexpr std::array<size_t, 4> make_qp_ind_array() {
  return {
      idx::qp_ind[0],
      idx::qp_ind[1],
      idx::qp_ind[2],
      idx::qp_ind[3]
  };
}

inline void precip(const real_t (&params)[3], real_t (&precip)[3], real_t zeta,
            real_t vc, real_t flx, real_t vt, real_t q, real_t q_kp1,
            real_t rho) {
  real_t rho_x, flx_eff, flx_partial;
  rho_x = q * rho;
  flx_eff = (rho_x / zeta) + static_cast<real_t>(2.0) * flx;
  flx_partial = rho_x * vc * fall_speed(rho_x, params);
  flx_partial = std::fmin(flx_partial, flx_eff);
  precip[0] = (zeta * (flx_eff - flx_partial)) /
              ((static_cast<real_t>(1.0) + zeta * vt) * rho); // q update
  precip[1] =
      (precip[0] * rho * vt + flx_partial) * static_cast<real_t>(0.5); // flx
  rho_x = (precip[0] + q_kp1) * static_cast<real_t>(0.5) * rho;
  precip[2] = vc * fall_speed(rho_x, params); // vt
}

void graupel(size_t &nvec, size_t &ke, size_t &ivstart, size_t &ivend,
             size_t &kstart, real_t &dt, array_1d_t<real_t> &dz,
             array_1d_t<real_t> &t, array_1d_t<real_t> &rho,
             array_1d_t<real_t> &p, array_1d_t<real_t> &qv,
             array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
             array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
             array_1d_t<real_t> &qg, real_t &qnc, array_1d_t<real_t> &prr_gsp,
             array_1d_t<real_t> &pri_gsp, array_1d_t<real_t> &prs_gsp,
             array_1d_t<real_t> &prg_gsp, array_1d_t<real_t> &pre_gsp,
             array_1d_t<real_t> &pflx) {

  // TODO: implement based on sequential code in implementations/sequential/graupel.cpp
  // std::cout << "sequential graupel" << std::endl;
  const auto qp_ind_copy = make_qp_ind_array();
  const auto qx_ind_copy = make_qx_ind_array(); 

  array_1d_t<size_t> kmin(nvec * np, ke + 1); // first level with condensate

  array_1d_t<real_t> emptyArray;
  // Arrays to store pointers for p and x
  array_1d_t<real_t*> p_array;
  array_1d_t<real_t*> x_array;

  // Add pointers to each array instead of storing the whole data
  p_array.emplace_back(prr_gsp.data());
  p_array.emplace_back(pri_gsp.data());
  p_array.emplace_back(prs_gsp.data());
  p_array.emplace_back(prg_gsp.data());
  p_array.emplace_back(emptyArray.data()); // Empty array pointer
  p_array.emplace_back(emptyArray.data());

  x_array.emplace_back(qr.data());
  x_array.emplace_back(qi.data());
  x_array.emplace_back(qs.data());
  x_array.emplace_back(qg.data());
  x_array.emplace_back(qc.data());
  x_array.emplace_back(qv.data());

  // The loop is intentionally i<nlev; since we are using an unsigned integer
  // data type, when i reaches 0, and you try to decrement further, (to -1), it
  // wraps to the maximum value representable by size_t.
  size_t total_size = ivend - ivstart;
  array_1d_t<size_t> indices_test(T_SIZE);

  std::iota(indices_test.begin(), indices_test.end(), 0);

  array_1d_t<size_t> flags(ke * (ivend - ivstart + 1), 0);
  array_1d_t<size_t> prefixsum(ke * (ivend - ivstart + 1));

  size_t jmx_ = std::transform_reduce(
    std::execution::par_unseq,
    indices_test.begin(), indices_test.end(),
    size_t(0), 
    std::plus<size_t>(),
    [=, kmin = kmin.data(), x_array = x_array.data(),
      t = t.data(), rho = rho.data(), qp_ind = qp_ind_copy,
      flags = flags.data()](size_t tid) {
      
      size_t count = 0;
      size_t update[4];
      size_t oned_vec_index;


      size_t base = total_size / T_SIZE;
      size_t rem  = total_size % T_SIZE;
      size_t chunk_size = base + (tid < rem) * 1;
      size_t start_cell = tid * base + std::min(tid, (size_t)rem);
      size_t j;
      for (size_t jj = 0; jj < chunk_size; ++jj) {
        j = tid + jj * T_SIZE;
        #pragma unroll np
        for (size_t ix = 0; ix < np; ++ix) {
          update[ix] = kmin[j * np + qp_ind[ix]];
        }
  
        for (size_t i = ke - 1; i < ke; --i) {        
          oned_vec_index = i * ivend + j;          
          if ((std::max({x_array[lqc][oned_vec_index], x_array[lqr][oned_vec_index], 
                x_array[lqs][oned_vec_index], x_array[lqi][oned_vec_index], 
                x_array[lqg][oned_vec_index]}) > qmin) || ((t[oned_vec_index] < tfrz_het2) && 
                (x_array[lqv][oned_vec_index] > qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index])))) {
                  flags[oned_vec_index] = 1;
                  count ++;
          }
          
          #pragma unroll np
          for (size_t ix = 0; ix < np; ++ix) {
            bool val = (x_array[qp_ind[ix]][oned_vec_index] > qmin);
            update[ix] =  val * i + (!val) * update[ix];
          }
        }
  
        #pragma unroll np
        for (size_t ix = 0; ix < np; ++ix) {
          kmin[j * np + qp_ind[ix]] = update[ix];
        }
      }
      return count;
    }
  );

  array_1d_t<size_t> ind_i(jmx_);
  array_1d_t<size_t> ind_j(jmx_);
  std::exclusive_scan(std::execution::par_unseq, flags.begin(), flags.end(), prefixsum.begin(), 0);

  std::for_each(std::execution::par_unseq, indices_test.begin(), indices_test.end(),
  [=, flags = flags.data(), prefixsum = prefixsum.data(), ind_i = ind_i.data(), ind_j = ind_j.data()](size_t tid) {
      size_t oned_vec_index;
      size_t base = total_size / T_SIZE;
      size_t rem  = total_size % T_SIZE;
      size_t chunk_size = base + (tid < rem) * 1;
      size_t start_cell = tid * base + std::min(tid, (size_t)rem);
      size_t j;
      for (size_t jj = 0; jj < chunk_size; ++jj) {
        j = tid + jj * T_SIZE;
        for (size_t i = ke - 1; i < ke; --i) {        
          oned_vec_index = i * ivend + j;
          if (flags[oned_vec_index]) {
            ind_j[prefixsum[oned_vec_index]] = j;
            ind_i[prefixsum[oned_vec_index]] = i;
          }
        }
      }
  });
  
  std::for_each(std::execution::par_unseq, indices_test.begin(), indices_test.end(),
    [=, p = p.data(), t = t.data(), x_array = x_array.data(), rho = rho.data(),
    ind_i = ind_i.data(), ind_j = ind_j.data(), qp_ind = qp_ind_copy, qx_ind = qx_ind_copy](size_t tid) {

        size_t base = jmx_ / T_SIZE;
        size_t rem  = jmx_ % T_SIZE;
        size_t chunk_size = base + (tid < rem) * 1;
        size_t start_cell = tid * base + std::min(tid, (size_t)rem);
        size_t j;
        for (size_t jj = 0; jj < chunk_size; ++jj) {
          j = tid + jj * T_SIZE;
            real_t cv, eta, qvsi, qice, qliq, qtot, dvsw, dvsw0, dvsi, n_ice,
            m_ice, x_ice, n_snow, l_snow, ice_dep, stot;
            real_t sx2x_sum;
            real_t sx2x[nx][nx] = {ZERO};
            real_t sink[nx],     // tendencies
            dqdt[nx];     // tendencies
            size_t k = ind_i[j];
            size_t iv = ind_j[j];
            size_t oned_vec_index = k * ivend + iv;
        
            bool is_sig_present = std::max({x_array[lqs][oned_vec_index], x_array[lqi][oned_vec_index], x_array[lqg][oned_vec_index]}) > qmin;
    
            dvsw = x_array[lqv][oned_vec_index] -
                  qsat_rho(t[oned_vec_index], rho[oned_vec_index]);
            qvsi = qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index]);
            dvsi = x_array[lqv][oned_vec_index] - qvsi;
            n_snow = snow_number(t[oned_vec_index], rho[oned_vec_index],
                                  x_array[lqs][oned_vec_index]);
            l_snow = snow_lambda(rho[oned_vec_index], x_array[lqs][oned_vec_index], n_snow);
        
            sx2x[lqc][lqr] = cloud_to_rain(t[oned_vec_index], x_array[lqc][oned_vec_index],
                                            x_array[lqr][oned_vec_index], qnc);
            sx2x[lqr][lqv] = rain_to_vapor(t[oned_vec_index], rho[oned_vec_index],
                                          x_array[lqc][oned_vec_index],
                                          x_array[lqr][oned_vec_index], dvsw, dt);
            sx2x[lqc][lqi] = cloud_x_ice(t[oned_vec_index], x_array[lqc][oned_vec_index],
                                          x_array[lqi][oned_vec_index], dt);
            sx2x[lqi][lqc] = -std::fmin(sx2x[lqc][lqi], ZERO);
            sx2x[lqc][lqi] = std::fmax(sx2x[lqc][lqi], ZERO);
            sx2x[lqc][lqs] = cloud_to_snow(t[oned_vec_index], x_array[lqc][oned_vec_index],
                                            x_array[lqs][oned_vec_index], n_snow, l_snow);
            sx2x[lqc][lqg] =
                cloud_to_graupel(t[oned_vec_index], rho[oned_vec_index],
                                x_array[lqc][oned_vec_index], x_array[lqg][oned_vec_index]);
        
            if (t[oned_vec_index] < tmelt) {
              n_ice = ice_number(t[oned_vec_index], rho[oned_vec_index]);
              m_ice = ice_mass(x_array[lqi][oned_vec_index], n_ice);
              x_ice = ice_sticking(t[oned_vec_index]);
        
              if (is_sig_present) {
                eta = deposition_factor(
                    t[oned_vec_index],
                    qvsi); // neglect cloud depth cor. from gcsp_graupel
                sx2x[lqv][lqi] = vapor_x_ice(x_array[lqi][oned_vec_index], m_ice, eta, dvsi,
                                            rho[oned_vec_index], dt);
                sx2x[lqi][lqv] = -std::fmin(sx2x[lqv][lqi], ZERO);
                sx2x[lqv][lqi] = std::fmax(sx2x[lqv][lqi], ZERO);
                ice_dep = std::fmin(sx2x[lqv][lqi], dvsi / dt);
        
                sx2x[lqi][lqs] = deposition_auto_conversion(x_array[lqi][oned_vec_index],
                                                            m_ice, ice_dep);
                sx2x[lqi][lqs] = sx2x[lqi][lqs] + ice_to_snow(x_array[lqi][oned_vec_index],
                                                              n_snow, l_snow, x_ice);
                sx2x[lqi][lqg] = ice_to_graupel(
                    rho[oned_vec_index], x_array[lqr][oned_vec_index],
                    x_array[lqg][oned_vec_index], x_array[lqi][oned_vec_index], x_ice);
                sx2x[lqs][lqg] =
                    snow_to_graupel(t[oned_vec_index], rho[oned_vec_index],
                                    x_array[lqc][oned_vec_index], x_array[lqs][oned_vec_index]);
                sx2x[lqr][lqg] = rain_to_graupel(
                    t[oned_vec_index], rho[oned_vec_index], x_array[lqc][oned_vec_index],
                    x_array[lqr][oned_vec_index], x_array[lqi][oned_vec_index],
                    x_array[lqs][oned_vec_index], m_ice, dvsw, dt);
              }
              sx2x[lqv][lqi] =
                  sx2x[lqv][lqi] +
                  ice_deposition_nucleation(t[oned_vec_index], x_array[lqc][oned_vec_index],
                                            x_array[lqi][oned_vec_index], n_ice, dvsi, dt);
            } else {
              sx2x[lqc][lqr] = sx2x[lqc][lqr] + sx2x[lqc][lqs] + sx2x[lqc][lqg];
              sx2x[lqc][lqs] = ZERO;
              sx2x[lqc][lqg] = ZERO;
              ice_dep = ZERO;
              eta = ZERO;
            }
        
            if (is_sig_present) {
              dvsw0 = x_array[lqv][oned_vec_index] - qsat_rho(tmelt, rho[oned_vec_index]);
              sx2x[lqv][lqs] =
                  vapor_x_snow(t[oned_vec_index], p[oned_vec_index],
                              rho[oned_vec_index], x_array[lqs][oned_vec_index], n_snow,
                              l_snow, eta, ice_dep, dvsw, dvsi, dvsw0, dt);
              sx2x[lqs][lqv] = -std::fmin(sx2x[lqv][lqs], ZERO);
              sx2x[lqv][lqs] = std::fmax(sx2x[lqv][lqs], ZERO);
              sx2x[lqv][lqg] = vapor_x_graupel(
                  t[oned_vec_index], p[oned_vec_index], rho[oned_vec_index],
                  x_array[lqg][oned_vec_index], dvsw, dvsi, dvsw0, dt);
              sx2x[lqg][lqv] = -std::fmin(sx2x[lqv][lqg], ZERO);
              sx2x[lqv][lqg] = std::fmax(sx2x[lqv][lqg], ZERO);
              sx2x[lqs][lqr] =
                  snow_to_rain(t[oned_vec_index], p[oned_vec_index],
                              rho[oned_vec_index], dvsw0, x_array[lqs][oned_vec_index]);
              sx2x[lqg][lqr] =
                  graupel_to_rain(t[oned_vec_index], p[oned_vec_index],
                                  rho[oned_vec_index], dvsw0, x_array[lqg][oned_vec_index]);
            }
        
            #pragma unroll nx
            for (size_t ix = 0; ix < nx; ix++) {
              sink[qx_ind[ix]] = ZERO;
              if ((is_sig_present) or (qx_ind[ix] == lqc) or (qx_ind[ix] == lqv) or
                  (qx_ind[ix] == lqr)) {
        
                #pragma unroll nx
                for (size_t i = 0; i < nx; i++) {
                  sink[qx_ind[ix]] = sink[qx_ind[ix]] + sx2x[qx_ind[ix]][i];
                }
                stot = x_array[qx_ind[ix]][oned_vec_index] / dt;
        
                if ((sink[qx_ind[ix]] > stot) &&
                    (x_array[qx_ind[ix]][oned_vec_index] > qmin)) {
                  real_t nextSink = ZERO;
        
                  #pragma unroll nx
                  for (size_t i = 0; i < nx; i++) {
                    sx2x[qx_ind[ix]][i] = sx2x[qx_ind[ix]][i] * stot / sink[qx_ind[ix]];
                    nextSink = nextSink + sx2x[qx_ind[ix]][i];
                  }
                  sink[qx_ind[ix]] = nextSink;
                }
              }
            }
            #pragma unroll nx
            for (size_t ix = 0; ix < nx; ix++) {
              sx2x_sum = 0;
              #pragma unroll nx
              for (size_t i = 0; i < nx; i++) {
                sx2x_sum = sx2x_sum + sx2x[i][qx_ind[ix]];
              }
              dqdt[qx_ind[ix]] = sx2x_sum - sink[qx_ind[ix]];
              x_array[qx_ind[ix]][oned_vec_index] = std::fmax(
                  ZERO, x_array[qx_ind[ix]][oned_vec_index] + dqdt[qx_ind[ix]] * dt);
            }
        
            qice = x_array[lqs][oned_vec_index] + x_array[lqi][oned_vec_index] + x_array[lqg][oned_vec_index];
            qliq = x_array[lqc][oned_vec_index] + x_array[lqr][oned_vec_index];
            qtot = x_array[lqv][oned_vec_index] + qice + qliq;
            cv = cvd + (cvv - cvd) * qtot + (clw - cvv) * qliq +
                (ci - cvv) * qice; // qtot? or qv?
            t[oned_vec_index] =
                t[oned_vec_index] +
                dt *
                    ((dqdt[lqc] + dqdt[lqr]) * (lvc - (clw - cvv) * t[oned_vec_index]) +
                    (dqdt[lqi] + dqdt[lqs] + dqdt[lqg]) *
                        (lsc - (ci - cvv) * t[oned_vec_index])) /
                    cv;
        }
  });

  size_t k_end = (lrain) ? ke : kstart - 1;

  std::for_each(std::execution::par_unseq, indices_test.begin(), indices_test.end(),
                [=, kmin = kmin.data(), x_array = x_array.data(), p_array = p_array.data(), dz = dz.data(),
                  pflx = pflx.data(), t = t.data(), pre_gsp = pre_gsp.data(),
                  rho = rho.data(), p = p.data(),
                  qp_ind = qp_ind_copy, qx_ind = qx_ind_copy] (size_t tid) {
      bool flag;

      size_t oned_vec_index, kp1;

      real_t vc, zeta, qice, qliq, e_int, xrho;
      real_t update[3];

      size_t base = total_size / T_SIZE;
      size_t rem  = total_size % T_SIZE;
      size_t chunk_size = base + (tid < rem) * 1;
      size_t start_cell = tid * base + std::min(tid, (size_t)rem);
      size_t iv;
      for (size_t jj = 0; jj < chunk_size; ++jj) {
        iv = tid + jj * T_SIZE;
        real_t vt[np] = {ZERO};
        real_t eflx = ZERO;
  
        const real_t params[4][3] = {
          {14.58, 0.111, 1.0e-12},
          {1.25, 0.160, 1.0e-12},
          {57.80, static_cast<real_t>(0.5) / static_cast<real_t>(3.0), 1.0e-12},
          {12.24, 0.217, 1.0e-08}};
  
        const size_t threshold = *std::min_element(kmin + (iv * np), kmin + (iv * np + np));
  
        for (size_t k = kstart; k < k_end; k++) {
  
          if (k < threshold)
            continue;
  
          oned_vec_index = k * ivend + iv;
  
          kp1 = std::min(ke - 1, k + 1);
  
          qliq = x_array[lqc][oned_vec_index] + x_array[lqr][oned_vec_index];
          qice = x_array[lqs][oned_vec_index] + x_array[lqi][oned_vec_index] + x_array[lqg][oned_vec_index];
  
          e_int =
              internal_energy(t[oned_vec_index], x_array[lqv][oned_vec_index], qliq,
                              qice, rho[oned_vec_index], dz[oned_vec_index]) +
              eflx;
          zeta = dt / (2.0 * dz[oned_vec_index]);
          xrho = std::sqrt(rho_00 / rho[oned_vec_index]);
          
          #pragma unroll np
          for (size_t ix = 0; ix < np; ix++) {
            if (k < kmin[iv * np + qp_ind[ix]])
              continue;
            vc = vel_scale_factor(qp_ind[ix], xrho, rho[oned_vec_index],
                                  t[oned_vec_index],
                                  x_array[qp_ind[ix]][oned_vec_index]);
            precip(params[qp_ind[ix]], update, zeta, vc, p_array[qp_ind[ix]][iv],
                    vt[ix], x_array[qp_ind[ix]][oned_vec_index],
                    x_array[qp_ind[ix]][kp1 * ivend + iv], rho[oned_vec_index]);
                    x_array[qp_ind[ix]][oned_vec_index] = update[0];
                    p_array[qp_ind[ix]][iv] = update[1];
            vt[ix] = update[2];
          }
  
          pflx[oned_vec_index] = p_array[lqs][iv] + p_array[lqi][iv] + p_array[lqg][iv];
          eflx =
              dt * (p_array[lqr][iv] * (clw * t[oned_vec_index] -
                                    cvd * t[kp1 * ivend + iv] - lvc) +
                    pflx[oned_vec_index] * (ci * t[oned_vec_index] -
                                            cvd * t[kp1 * ivend + iv] - lsc));
          pflx[oned_vec_index] = pflx[oned_vec_index] + p_array[lqr][iv];
          qliq = x_array[lqc][oned_vec_index] + x_array[lqr][oned_vec_index];
          qice = x_array[lqs][oned_vec_index] + x_array[lqi][oned_vec_index] +
          x_array[lqg][oned_vec_index];
          e_int = e_int - eflx;
          t[oned_vec_index] =
              T_from_internal_energy(e_int, x_array[lqv][oned_vec_index], qliq, qice,
                                      rho[oned_vec_index], dz[oned_vec_index]);
          flag = (k == ke - 1);
          pre_gsp[iv] = (flag) * eflx / dt + (!flag) * pre_gsp[iv];
        }
      }
  });
}