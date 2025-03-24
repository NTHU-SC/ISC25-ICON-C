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
}
