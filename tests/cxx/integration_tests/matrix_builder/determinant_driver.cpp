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
using pt = simde::eval_braket<WFType, simde::type::op_base_type, WFType>;

using simde::type::tensor;

template<typename WFType>
using erased_type =
  chemist::braket::BraKet<WFType, simde::type::op_base_type, WFType>;

TEST_CASE("DeterminantDriver") {
    auto mm  = test_scf::load_modules();
    auto mod = mm.at("Determinant driver");

    using wf_type   = simde::type::rscf_wf;
    using index_set = typename wf_type::orbital_index_set_type;

    wf_type psi(index_set{0}, test_scf::h2_cmos());
    simde::type::many_electrons es{2};

    SECTION("Calling Kinetic") {
        simde::type::T_e_type T_e(es);
        chemist::braket::BraKet braket(psi, T_e, psi);
        erased_type<wf_type> copy_braket(braket);
        const auto& T = mod.run_as<pt<wf_type>>(copy_braket);
        REQUIRE_THAT(T, WithinAbs(0.637692, 1E-6));
    }

    SECTION("Calling Electron-Nuclear Attraction") {
        auto h2_nuclei = test_scf::make_h2<simde::type::nuclei>();
        simde::type::V_en_type V_en(es, h2_nuclei);
        chemist::braket::BraKet braket(psi, V_en, psi);
        erased_type<wf_type> copy_braket(braket);
        const auto& V = mod.run_as<pt<wf_type>>(copy_braket);
        REQUIRE_THAT(V, WithinAbs(-1.96830777255516853, 1E-6));
    }

    SECTION("Calling J") {
        auto rho = test_scf::h2_density();
        simde::type::J_e_type J_e(es, rho);
        chemist::braket::BraKet braket(psi, J_e, psi);
        erased_type<wf_type> copy_braket(braket);
        const auto& J = mod.run_as<pt<wf_type>>(copy_braket);
        REQUIRE_THAT(J, WithinAbs(0.76056339681664897, 1E-6));
    }

    SECTION("Calling K") {
        auto rho = test_scf::h2_density();
        simde::type::K_e_type K_e(es, rho);
        chemist::braket::BraKet braket(psi, K_e, psi);
        erased_type<wf_type> copy_braket(braket);
        const auto& K = mod.run_as<pt<wf_type>>(copy_braket);
        REQUIRE_THAT(K, WithinAbs(-0.76056339681664897, 1E-6));
    }
}