/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "gauxc.hpp"

namespace scf::xc::gauxc {
namespace {

const auto desc = R"(
Chemist AOBasisSet to GauXC BasisSet
====================================

TODO: Write me!!
)";

}

// Module to convert NWX Basis -> GauXC Basis
MODULE_CTOR(BasisConversion) {
    satisfies_property_type<gauxc_basis_conversion_t>();
    description(desc);
}
MODULE_RUN(BasisConversion) {
    const auto& [nwx_basis] = gauxc_basis_conversion_t::unwrap_inputs(inputs);

    GauXC::BasisSet<double> gauxc_basis;
    auto nshells = nwx_basis.n_shells();
    gauxc_basis.reserve(nshells);
    for(decltype(nshells) i_sh = 0; i_sh < nshells; ++i_sh) {
        const auto shell      = nwx_basis.shell(i_sh);
        const auto nprims     = shell.n_primitives();
        const auto first_prim = shell.primitive(0);
        const auto last_prim  = shell.primitive(nprims - 1);
        const int l           = shell.l();
        const bool pure       = shell.pure() == chemist::ShellType::pure;

        if(nprims > GauXC::detail::shell_nprim_max) {
            std::stringstream ss;
            ss << "GauXC received a shell with nprim = " << nprims
               << "with NPRIM_MAX = " << GauXC::detail::shell_nprim_max
               << ". Please reconfigure GauXC.";
            throw std::runtime_error(ss.str());
        }

        GauXC::Shell<double>::prim_array alpha, coeff;
        GauXC::Shell<double>::cart_array center = {
          shell.center().x(), shell.center().y(), shell.center().z()};

        std::copy(&first_prim.exponent(), &last_prim.exponent() + 1,
                  alpha.begin());
        std::copy(&first_prim.coefficient(), &last_prim.coefficient() + 1,
                  coeff.begin());

        gauxc_basis.emplace_back(
          GauXC::PrimSize(nprims), GauXC::AngularMomentum(l),
          GauXC::SphericalType(pure), alpha, coeff, center);
    }

    auto rv = results();
    return gauxc_basis_conversion_t::wrap_results(rv, gauxc_basis);
}
} // namespace scf::xc::gauxc
