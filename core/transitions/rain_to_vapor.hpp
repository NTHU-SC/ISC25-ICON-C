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

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace transition {

/**
 * @brief TODO
 *
 * @param [in] t Temperature
 * @param [in] rho Ambient density
 * @param [in] qc Specific humidity of cloud
 * @param [in] qr Specific humidity of rain
 * @param [in] dvsw qv - qsat_water(T)
 * @param [in] dt Time step
 * @return Mass from qc to qr
 */
TARGET real_t rain_to_vapor(real_t t, real_t rho, real_t qc, real_t qr,
                            real_t dvsw, real_t dt) {
  constexpr real_t b1_rv =
      real_t{0.16667}; // exponent in power-law relation for mass density
  constexpr real_t b2_rv = real_t{0.55555};
  constexpr real_t c1_rv = real_t{0.61};
  constexpr real_t c2_rv = real_t{-0.0163};
  constexpr real_t c3_rv = real_t{1.111e-4};
  constexpr real_t a1_rv = real_t{1.536e-3};
  constexpr real_t a2_rv = real_t{1.0E+0}; // constant in rain evap formula
  constexpr real_t a3_rv = real_t{
      19.0621E+0}; // prefactor (from gamma dist. and properties of air/water)
  if (qr > graupel_ct::qmin && (dvsw + qc <= real_t{0.0})) {
    real_t tc = t - thermodyn::tmelt;
    real_t evap_max = (c1_rv + tc * (c2_rv + c3_rv * tc)) * (-dvsw) / dt;
    return fmin(a1_rv * (a2_rv + a3_rv * pow((qr * rho), b1_rv)) * (-dvsw) *
                    pow((qr * rho), b2_rv),
                evap_max);
  }

  return static_cast<real_t>(0.0);
}

} // namespace transition
