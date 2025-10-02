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

using grid_pt = simde::MolecularGrid;
using ao_pt   = simde::AOCollocationMatrix;
using rho_pt  = simde::EDensityCollocationMatrix;
using pt      = simde::aos_xc_e_aos;

TEST_CASE("LibXCPotential") {
    using float_type = double;
    using tensorwrapper::operations::approximately_equal;

    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    mm.change_submod("LibXC Potential", "Integration grid", "Grid From File");
    auto& mod = mm.at("LibXC Potential");

    auto rho  = test_scf::h2_density<float_type>();
    auto aos  = test_scf::h2_aos();
    auto func = chemist::qm_operator::xc_functional::SVWN3;

    // For this test only the weights matter (we hardcode the density array)
    using grid_point = chemist::GridPoint;
    std::vector<grid_point> grid_points;
    grid_points.push_back(grid_point(15.80337241, 0.0, 0.0, 0.0));
    grid_points.push_back(grid_point(0.03961399, 0.0, 0.0, 0.0));
    grid_points.push_back(grid_point(15.80337241, 0.0, 0.0, -1.3984));
    grid_points.push_back(grid_point(0.03961399, 0.0, 0.0, -1.3984));
    chemist::Grid grid(grid_points.begin(), grid_points.end());

    auto grid_mod = pluginplay::make_lambda<grid_pt>([&](auto&& mol) {
        REQUIRE(mol.size() == 2);
        return grid;
    });

    auto rho_mod =
      pluginplay::make_lambda<rho_pt>([&](auto&& grid_in, auto&& P) {
          REQUIRE(grid_in == grid);
          REQUIRE(P == rho);
          return simde::type::tensor(
            {0.01313525, 0.32141338, 0.01313525, 0.3214133});
      });

    mod.change_submod("Integration grid", grid_mod);
    mod.change_submod("Density on a grid", rho_mod);

    simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
    chemist::braket::BraKet braket(aos, xc_op, aos);

    auto vxc = mod.run_as<pt>(braket);
    simde::type::tensor corr{{-1.83222, -0.402489}, {-0.402489, -0.0886007}};
    REQUIRE(approximately_equal(vxc, corr, 1e-5));
}
