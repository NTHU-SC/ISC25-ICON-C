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
 * @brief Conversion rate of ice to snow
 *
 * @param [in] qi Ice specific mass
 * @param [in] ns Snow number
 * @param [in] lambda Snow intercept parameter, lambda
 * @param [in] sticking_eff Ice sticking effiency
 * @return conversion rate of ice to snow
 */
TARGET real_t ice_to_snow(real_t qi, real_t ns, real_t lambda,
                          real_t sticking_eff) {

  constexpr real_t qi0 =
      real_t{0.0}; // critical ice required for autoconversion
  constexpr real_t c_iau = real_t{1.0E-3}; // coefficient of auto conversion
  constexpr real_t c_agg =
      static_cast<real_t>(2.61) *
      graupel_ct::v0s; // coeff of aggregation (2.610 = pi*gam(v1s+3)/4)
  constexpr real_t b_agg =
      -(graupel_ct::v1s + static_cast<real_t>(3.0)); // aggregation exponent

  return (qi > graupel_ct::qmin)
             ? sticking_eff *
                   (c_iau * fmax(static_cast<real_t>(0.0), (qi - qi0)) +
                    qi * (c_agg * ns) * pow(lambda, b_agg))
             : static_cast<real_t>(0.0);
}
} // namespace transition
