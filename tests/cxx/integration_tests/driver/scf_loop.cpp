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

TEMPLATE_LIST_TEST_CASE("SCFLoop", "", test_scf::float_types) {
    using float_type = TestType;
    using wf_type    = simde::type::rscf_wf;
    using index_set  = typename wf_type::orbital_index_set_type;

    auto mm   = test_scf::load_modules<float_type>();
    auto& mod = mm.at("Loop");

    using tensorwrapper::buffer::make_contiguous;
    tensorwrapper::shape::Smooth shape_corr{};
    auto pcorr = make_contiguous<float_type>(shape_corr);
    using tensorwrapper::operations::approximately_equal;

    SECTION("H2") {
        wf_type psi0(index_set{0}, test_scf::h2_cmos<float_type>());

        auto H = test_scf::h2_hamiltonian();
        chemist::braket::BraKet H_00(psi0, H, psi0);

        SECTION("Default") {
            mod.change_input("DIIS", true);

            const auto& [e, psi] = mod.template run_as<pt<wf_type>>(H_00, psi0);
            pcorr.set_elem({}, -1.1167592336);
            tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
            REQUIRE(approximately_equal(corr, e, 1E-6));
        }

        SECTION("With DIIS") {
            const auto& [e, psi] = mod.template run_as<pt<wf_type>>(H_00, psi0);
            pcorr.set_elem({}, -1.1167592336);
            tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
            REQUIRE(approximately_equal(corr, e, 1E-6));
        }
    }

    SECTION("He") {
        wf_type psi0(index_set{0}, test_scf::he_cmos<float_type>());

        auto H = test_scf::he_hamiltonian();
        chemist::braket::BraKet H_00(psi0, H, psi0);

        SECTION("Default") {
            mod.change_input("DIIS", true);
            const auto& [e, psi] = mod.template run_as<pt<wf_type>>(H_00, psi0);

            pcorr.set_elem({}, -2.807783957539);
            tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
            REQUIRE(approximately_equal(corr, e, 1E-6));
        }

        SECTION("With DIIS") {
            const auto& [e, psi] = mod.template run_as<pt<wf_type>>(H_00, psi0);

            pcorr.set_elem({}, -2.807783957539);
            tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
            REQUIRE(approximately_equal(corr, e, 1E-6));
        }
    }
}
