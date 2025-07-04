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
 * @param [in] qc Snow specific mass
 * @param [in] qg Graupel specific mass
 * @return Graupel riming rate
 */
TARGET real_t cloud_to_graupel(real_t t, real_t rho, real_t qc, real_t qg) {

  constexpr real_t a_rim = real_t{4.43};
  constexpr real_t b_rim = real_t{0.94878};
  return (fmin(qc, qg) > graupel_ct::qmin && t > graupel_ct::tfrz_hom)
             ? a_rim * qc * pow(qg * rho, b_rim)
             : static_cast<real_t>(0.0);
}
} // namespace transition
