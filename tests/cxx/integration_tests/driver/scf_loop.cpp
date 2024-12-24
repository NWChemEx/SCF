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

#include "../integration_tests.hpp"

using Catch::Matchers::WithinAbs;

template<typename WFType>
using egy_pt = simde::eval_braket<WFType, simde::type::hamiltonian, WFType>;

template<typename WFType>
using pt = simde::Optimize<egy_pt<WFType>, WFType>;

TEST_CASE("SCFLoop") {
    using wf_type = simde::type::rscf_wf;
    auto mm       = test_scf::load_modules();
    auto& mod     = mm.at("Loop");

    using index_set = typename wf_type::orbital_index_set_type;
    wf_type psi0(index_set{0}, test_scf::h2_cmos());

    auto H = test_scf::h2_hamiltonian();

    chemist::braket::BraKet H_00(psi0, H, psi0);
    const auto& [e, psi] = mod.run_as<pt<wf_type>>(H_00, psi0);
    REQUIRE_THAT(e, WithinAbs(-1.1167592336, 1E-6));
}