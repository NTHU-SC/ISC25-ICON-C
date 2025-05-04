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
 * @param [in] qi Specific humidity of ice
 * @param [in] mi Ice crystal mass
 * @param [in] eta Deposition factor
 * @param [in] dvsi Vapor excess with respect to ice sat
 * @param [in] dt Time step
 * @return Rate of vapor deposition to ice
 */
TARGET real_t vapor_x_ice(real_t qi, real_t mi, real_t eta, real_t dvsi,
                          real_t rho, real_t dt) {
  constexpr real_t ami =
      real_t{130.0}; // Formfactor for mass-size relation of cld ice
  constexpr real_t b_exp =
      real_t{-0.67}; // exp. for conv. (-1 + 0.33) of ice mass to sfc area
  const real_t a_fact =
      static_cast<real_t>(4.0) *
      pow(ami, static_cast<real_t>(-1.0) / static_cast<real_t>(3.0));
  real_t result = real_t{0.0};

  if (qi > graupel_ct::qmin) {
    result = (a_fact * eta) * rho * qi * pow(mi, b_exp) * dvsi;

    if (result > 0.) {
      result = fmin(result, dvsi / dt);
    } else {
      result = fmax(result, dvsi / dt);
      result = fmax(result, -qi / dt);
    }
  }

  return result;
}
} // namespace transition
