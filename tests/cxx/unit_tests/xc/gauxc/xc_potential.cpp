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

using pt = simde::aos_xc_e_aos;
using tensorwrapper::operations::approximately_equal;

TEST_CASE("XCPotential") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod         = mm.at("GauXC XC Potential");
    using nuclei_type = simde::type::nuclei;

    simde::type::electron e;
    auto func = chemist::qm_operator::xc_functional::PBE0;

    SECTION("h2") {
        auto rho        = test_scf::h2_density<double>();
        const auto& aos = rho.basis_set();
        simde::type::xc_e_type xc_op(func, e, rho);
        simde::type::braket xc_ij(aos, xc_op, aos);
        auto vxc = mod.run_as<pt>(xc_ij);

        simde::type::tensor corr{{-0.357302, -0.23347}, {-0.23347, -0.357302}};
        REQUIRE(approximately_equal(vxc, corr, 1E-5));
    }

    SECTION("he") {
        auto rho        = test_scf::he_density<double>();
        const auto& aos = rho.basis_set();
        simde::type::xc_e_type xc_op(func, e, rho);
        simde::type::braket xc_ij(aos, xc_op, aos);
        auto vxc = mod.run_as<pt>(xc_ij);

        tensorwrapper::allocator::Eigen<double> alloc(mm.get_runtime());
        tensorwrapper::shape::Smooth shape_corr{1, 1};
        auto pcorr =
          alloc.allocate(tensorwrapper::layout::Physical(shape_corr));
        pcorr->set_elem({0, 0}, -0.526535);
        simde::type::tensor corr(shape_corr, std::move(pcorr));
        REQUIRE(approximately_equal(vxc, corr, 1E-5));
    }
}