/*
 * Copyright 2026 NWChemEx-Project
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

using grid_pt = simde::MolecularGrid;
using pt      = simde::aos_xc_e_aos;
using quad_pt = scf::property_types::xc_potential_quadrature;

// The wrapper should build the grid and forward it, with its inputs, to the
// quadrature submodule, returning the quadrature's result unchanged.
TEST_CASE("XCPotentialViaQuadrature") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();

    auto rho  = test_scf::h2_density<double>();
    auto aos  = test_scf::h2_aos();
    auto func = chemist::qm_operator::xc_functional::SVWN3;
    simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
    chemist::braket::BraKet braket(aos, xc_op, aos);
    auto grid = test_scf::h2_xc_grid();
    simde::type::tensor corr{{1.0, 2.0}, {3.0, 4.0}};

    auto grid_mod = pluginplay::make_lambda<grid_pt>([&](auto&& mol) {
        REQUIRE(mol.size() == 2);
        return grid;
    });

    auto quad_mod =
      pluginplay::make_lambda<quad_pt>([&](auto&& braket_in, auto&& grid_in) {
          REQUIRE(braket_in == braket);
          REQUIRE(grid_in == grid);
          return corr;
      });

    for(std::string key : {"SVWN3 Potential", "LibXC Potential"}) {
        auto& mod = mm.at(key);
        mod.change_submod("Integration grid", grid_mod);
        mod.change_submod("Quadrature", quad_mod);
        REQUIRE(mod.run_as<pt>(braket) == corr);
    }
}
