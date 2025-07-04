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

namespace property {

/**
 * @brief Ice number following cooper
 *
 * @param [in] t Ambient temperature (kelvin)
 * @param [in] rho Ambient density
 * @return Ice number
 */
TARGET real_t ice_number(real_t t, real_t rho) {
  constexpr real_t a_coop = real_t{5.000};  // parameter in cooper fit
  constexpr real_t b_coop = real_t{0.304};  // parameter in cooper fit
  constexpr real_t nimax = real_t{250.e+3}; // maximal number of ice crystals
  return fmin(nimax, a_coop * exp(b_coop * (thermodyn::tmelt - t))) / rho;
}

} // namespace property
