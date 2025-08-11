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

#include "../../../test_scf.hpp"
#include <iostream>
#include <pluginplay/pluginplay.hpp>

using namespace scf;

using XC_e_t = simde::type::XC_e_type;

template<typename WFType>
using pt = simde::eval_braket<WFType, XC_e_t, WFType>;

using tensorwrapper::operations::approximately_equal;

TEST_CASE("XCEnergy") {
    using wf_type    = simde::type::rscf_wf;
    using index_set  = typename wf_type::orbital_index_set_type;
    using float_type = double;

    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("GauXC XC Energy");

    auto func = chemist::qm_operator::xc_functional::PBE0;

    SECTION("h2") {
        wf_type psi(index_set{0}, test_scf::h2_cmos<float_type>());
        auto rho = test_scf::h2_density<float_type>();
        simde::type::many_electrons es{2};
        XC_e_t xc_op(func, es, rho);
        chemist::braket::BraKet xc_ij(psi, xc_op, psi);
        auto exc = mod.run_as<pt<wf_type>>(xc_ij);
        simde::type::tensor corr(-0.587164);
        REQUIRE(approximately_equal(exc, corr, 1E-5));
    }

    SECTION("he") {
        wf_type psi(index_set{0}, test_scf::he_cmos<float_type>());
        auto rho = test_scf::he_density<double>();
        simde::type::many_electrons es{2};
        XC_e_t xc_op(func, es, rho);
        chemist::braket::BraKet xc_ij(psi, xc_op, psi);
        auto exc = mod.run_as<pt<wf_type>>(xc_ij);
        simde::type::tensor corr(-0.819986);
        REQUIRE(approximately_equal(exc, corr, 1E-5));
    }
}
