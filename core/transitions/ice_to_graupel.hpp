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
 * @param [in] rho Density
 * @param [in] qr Rain specific mass
 * @param [in] qg Graupel specific mass
 * @param [in] qi Ice specific mass
 * @param [in] sticking_eff Sticking effiency
 * @return aggregation of ice by graupel
 */
TARGET real_t ice_to_graupel(real_t rho, real_t qr, real_t qg, real_t qi,
                             real_t sticking_eff) {
  constexpr real_t a_ct =
      real_t{1.72}; //  (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
  constexpr real_t b_ct = static_cast<real_t>(7.0) / static_cast<real_t>(8.0);
  constexpr real_t c_agg_ct = real_t{2.46};
  constexpr real_t b_agg_ct = real_t{0.94878};
  real_t result = real_t{0.0};
  if (qi > graupel_ct::qmin) {
    if (qg > graupel_ct::qmin) {
      result = sticking_eff * qi * c_agg_ct * pow(rho * qg, b_agg_ct);
    }
    if (qr > graupel_ct::qmin) {
      result = result + a_ct * qi * pow(rho * qr, b_ct);
    }
  }

  return result;
}

} // namespace transition
