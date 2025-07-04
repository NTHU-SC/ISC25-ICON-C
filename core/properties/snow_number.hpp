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

constexpr real_t n0s0 = real_t{8.00e+5};
/**
 * @brief TODO
 * @param [in] t Temperature
 * @param [in] rho Ambient air density
 * @param [in] qs Snow specific mass
 * @return Snow number
 */
TARGET real_t snow_number(real_t t, real_t rho, real_t qs) {

  constexpr real_t tmin = thermodyn::tmelt - static_cast<real_t>(40.);
  constexpr real_t tmax = thermodyn::tmelt;
  constexpr real_t qsmin = real_t{2.0e-6};
  constexpr real_t xa1 = real_t{-1.65e+0};
  constexpr real_t xa2 = real_t{5.45e-2};
  constexpr real_t xa3 = real_t{3.27e-4};
  constexpr real_t xb1 = real_t{1.42e+0};
  constexpr real_t xb2 = real_t{1.19e-2};
  constexpr real_t xb3 = real_t{9.60e-5};
  constexpr real_t n0s1 =
      static_cast<real_t>(13.5) * static_cast<real_t>(5.65e+05);
  constexpr real_t n0s2 = real_t{-0.107};
  constexpr real_t n0s3 = real_t{13.5};
  constexpr real_t n0s4 = static_cast<real_t>(0.5) * n0s1;
  constexpr real_t n0s5 = real_t{1.e6};
  constexpr real_t n0s6 = static_cast<real_t>(1.e2) * n0s1;
  constexpr real_t n0s7 = real_t{1.e9};

  if (qs > graupel_ct::qmin) {
    real_t tc = fmax(fmin(t, tmax), tmin) - thermodyn::tmelt;
    real_t alf = pow(static_cast<real_t>(10.), (xa1 + tc * (xa2 + tc * xa3)));
    real_t bet = xb1 + tc * (xb2 + tc * xb3);
    real_t n0s =
        n0s3 *
        pow(((qs + qsmin) * rho / graupel_ct::ams),
            (static_cast<real_t>(4.0) - static_cast<real_t>(3) * bet)) /
        (alf * alf * alf);
    real_t y = exp(n0s2 * tc);
    real_t n0smn = fmax(n0s4 * y, n0s5);
    real_t n0smx = fmin(n0s6 * y, n0s7);
    return fmin(n0smx, fmax(n0smn, n0s));
  } else {
    return n0s0;
  }
}
} // namespace property
