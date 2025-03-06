/*
 * Copyright 2024 NWChemEx-Project
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

#include "matrix_builder.hpp"

namespace scf::matrix_builder {

namespace {

const auto desc = R"(
)";

}

using simde::type::electronic_hamiltonian;

template<typename WFType>
using pt = simde::eval_braket<WFType, electronic_hamiltonian, WFType>;

template<typename WFType>
using det_pt = simde::eval_braket<WFType, simde::type::op_base_type, WFType>;

MODULE_CTOR(ElectronicEnergy) {
    using wf_type = simde::type::rscf_wf;
    description(desc);
    satisfies_property_type<pt<wf_type>>();
    add_submodule<det_pt<wf_type>>("determinant driver");
}

MODULE_RUN(ElectronicEnergy) {
    using wf_type       = simde::type::rscf_wf;
    const auto&& [H_ij] = pt<wf_type>::unwrap_inputs(inputs);
    const auto& bra     = H_ij.bra();
    const auto& H       = H_ij.op();
    const auto& ket     = H_ij.ket();

    auto& mod = submods.at("determinant driver");
    simde::type::tensor e;

    auto n_ops = H.size();
    for(decltype(n_ops) i = 0; i < n_ops; ++i) {
        const auto& ci  = H.coefficient(i);
        const auto& O_i = H.get_operator(i);

        chemist::braket::BraKet O_ij(bra, O_i, ket);
        const auto o = mod.run_as<det_pt<wf_type>>(O_ij);
        simde::type::tensor temp;
        temp("") = o("") * ci;
        if(i == 0)
            e = temp;
        else { e("") = e("") + temp(""); }
    }

    auto rv = results();
    return pt<wf_type>::wrap_results(rv, e);
}

} // namespace scf::matrix_builder