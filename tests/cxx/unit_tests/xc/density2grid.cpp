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

#include "../../test_scf.hpp"
#include <iostream>
#include <pluginplay/pluginplay.hpp>

using namespace scf;

using pt    = simde::EDensityCollocationMatrix;
using ao_pt = simde::AOCollocationMatrix;
using tensorwrapper::operations::approximately_equal;

TEST_CASE("Density2Grid") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("Density2Grid");

    chemist::Grid grid; // Value doesn't matter for this test b/c of submodule

    SECTION("H2") {
        auto rho = test_scf::h2_density<double>();
        simde::type::tensor ao2grid{{1.0, 2.0, 3.0, 4.0}, {5.0, 6.0, 7.0, 8.0}};
        auto ao_mod =
          pluginplay::make_lambda<ao_pt>([&](auto&& grid_in, auto&& aos) {
              REQUIRE(grid_in == grid);
              REQUIRE(aos == rho.basis_set().ao_basis_set());
              return ao2grid;
          });
        mod.change_submod("AOs on a grid", ao_mod);
        auto rv = mod.run_as<pt>(grid, rho);
        simde::type::tensor corr{11.513088, 20.467712, 31.9808, 46.052352};
        REQUIRE(approximately_equal(rv, corr, 1e-4));
    }
}
