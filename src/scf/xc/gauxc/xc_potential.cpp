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
Exchange-Correlation Potential via Driver
=========================================
)";
}

using pt = simde::aos_xc_e_aos;

MODULE_CTOR(XCPotential) {
    satisfies_property_type<pt>();
    description(desc);

    add_submodule<XCDriver>("XC Driver");
}

MODULE_RUN(XCPotential) {
    const auto& [braket] = pt::unwrap_inputs(inputs);

    const auto& bra_aos = braket.bra();
    const auto& xc_op   = braket.op();
    const auto& ket_aos = braket.ket();

    const auto func = xc_op.functional_name();
    const auto& P   = xc_op.rhs_particle();
    const auto& aos = bra_aos.ao_basis_set();

    auto& driver           = submods.at("XC Driver");
    const auto& [exc, vxc] = driver.run_as<XCDriver>(func, aos, P.value());

    auto rv = results();
    return pt::wrap_results(rv, vxc);
}

} // namespace scf::xc::gauxc