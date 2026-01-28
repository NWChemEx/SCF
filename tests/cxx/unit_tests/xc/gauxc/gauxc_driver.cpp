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

#include "../../../test_scf.hpp"
#include <iostream>
#include <pluginplay/pluginplay.hpp>
#include <scf/xc/gauxc/gauxc_property_types.hpp>

using namespace scf;

using pt = scf::xc::gauxc::XCDriver;
using tensorwrapper::operations::approximately_equal;

TEST_CASE("GauXCDriver") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("GauXC Driver");

    auto func = chemist::qm_operator::xc_functional::PBE0;

    SECTION("h2") {
        auto rho        = test_scf::h2_density<double>();
        const auto& aos = rho.basis_set().ao_basis_set();
        auto [exc, vxc] = mod.run_as<pt>(func, aos, rho.value());
        simde::type::tensor exc_corr(-0.587164);
        REQUIRE(approximately_equal(exc_corr, exc, 1E-5));

        simde::type::tensor vxc_corr{{-0.357302, -0.23347},
                                     {-0.23347, -0.357302}};
        REQUIRE(approximately_equal(vxc, vxc_corr, 1E-5));
    }

    SECTION("he") {
        auto rho        = test_scf::he_density<double>();
        const auto& aos = rho.basis_set().ao_basis_set();

        auto [exc, vxc] = mod.run_as<pt>(func, aos, rho.value());
        simde::type::tensor exc_corr(-0.819986);
        REQUIRE(approximately_equal(exc_corr, exc, 1E-5));

        using tensorwrapper::buffer::make_contiguous;
        tensorwrapper::shape::Smooth shape_corr{1, 1};
        auto pcorr = make_contiguous<double>(shape_corr);
        pcorr.set_elem({0, 0}, -0.526535);
        simde::type::tensor vxc_corr(shape_corr, std::move(pcorr));
        REQUIRE(approximately_equal(vxc, vxc_corr, 1E-5));
    }
}
