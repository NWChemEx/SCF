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

#include "../test_scf.hpp"
#include <scf/scf.hpp>
#include <scf/xc/gauxc/gauxc_property_types.hpp>

using pt = scf::xc::gauxc::gauxc_basis_conversion_t;

TEST_CASE("BasisConversion") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("GauXC Basis Converter");

    using nuclei_type = simde::type::nuclei;

    SECTION("he") {
        auto mol = test_scf::make_he<nuclei_type>();
        auto bs  = test_scf::he_basis(mol);

        auto rv = mod.run_as<pt>(bs);
        REQUIRE(rv.nbf() == 1);
        REQUIRE(rv.nshells() == 1);
        REQUIRE(rv.max_l() == 0);
    }

    SECTION("h2") {
        auto mol = test_scf::make_h2<nuclei_type>();
        auto bs  = test_scf::h_basis(mol);

        auto rv = mod.run_as<pt>(bs);
        REQUIRE(rv.nbf() == 2);
        REQUIRE(rv.nshells() == 2);
        REQUIRE(rv.max_l() == 0);
    }
}