/*
 * Copyright 2025 NWChemEx-Project
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

using XC_e_t = simde::type::XC_e_type;

template<typename WFType>
using pt = simde::eval_braket<WFType, XC_e_t, WFType>;

MODULE_CTOR(XCEnergy) {
    using wf_type = simde::type::rscf_wf;
    satisfies_property_type<pt<wf_type>>();

    add_submodule<XCDriver>("XC Driver");
}

MODULE_RUN(XCEnergy) {
    using wf_type        = simde::type::rscf_wf;
    const auto& [braket] = pt<wf_type>::unwrap_inputs(inputs);

    const auto& bra_wf = braket.bra();
    const auto& xc_op  = braket.op();
    const auto& ket_wf = braket.ket();

    if(bra_wf != ket_wf)
        throw std::runtime_error("Expected the same basis set");

    const auto func = xc_op.functional_name();
    const auto& P   = xc_op.rhs_particle();
    const auto& aos = P.basis_set().ao_basis_set();

    auto& driver           = submods.at("XC Driver");
    const auto& [exc, vxc] = driver.run_as<XCDriver>(func, aos, P.value());

    auto rv = results();
    return pt<wf_type>::wrap_results(rv, exc);
}

} // namespace scf::xc::gauxc
