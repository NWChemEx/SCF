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

using grid_pt  = simde::MolecularGrid;
using rho_pt   = simde::EDensityCollocationMatrix;
using pt       = simde::aos_xc_e_aos;
using quad_pt  = scf::property_types::xc_potential_quadrature;
using error_pt = scf::property_types::xc_potential_uncertainty;
using tensorwrapper::operations::approximately_equal;

TEST_CASE("UQXCPotential") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("UQ XC Potential");

    auto rho  = test_scf::h2_density<double>();
    auto aos  = test_scf::h2_aos();
    auto func = chemist::qm_operator::xc_functional::SVWN3;
    simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
    chemist::braket::BraKet braket(aos, xc_op, aos);
    auto grid = test_scf::h2_xc_grid();

    auto grid_mod = pluginplay::make_lambda<grid_pt>([&](auto&& mol) {
        REQUIRE(mol.size() == 2);
        return grid;
    });
    mod.change_submod("Integration grid", grid_mod);

    // Stands in for the uncertain result the error model returns
    simde::type::tensor uq_corr{{5.0, 6.0}, {7.0, 8.0}};

    SECTION("Forwards the grid and the quadrature to the error model") {
        simde::type::tensor vxc{{1.0, 2.0}, {3.0, 4.0}};
        auto quad_mod = pluginplay::make_lambda<quad_pt>(
          [&](auto&& braket_in, auto&& grid_in) {
              REQUIRE(braket_in == braket);
              REQUIRE(grid_in == grid);
              return vxc;
          });

        auto error_mod = pluginplay::make_lambda<error_pt>(
          [&](auto&& braket_in, auto&& grid_in, auto&& vxc_in) {
              REQUIRE(braket_in == braket);
              REQUIRE(grid_in == grid);
              REQUIRE(vxc_in == vxc);
              return uq_corr;
          });

        mod.change_submod("XC Potential", quad_mod);
        mod.change_submod("Incompleteness Error Model", error_mod);
        REQUIRE(mod.run_as<pt>(braket) == uq_corr);
    }

    SECTION("Defaults to the native SVWN3 quadrature") {
        // Hardcode the density on the grid (see the SVWN3PotentialQuadrature
        // test, whose reference values these are)
        auto rho_mod = pluginplay::make_lambda<rho_pt>(
          [&](auto&&, auto&&) { return test_scf::h2_rho_on_xc_grid(); });
        mm.at("SVWN3 Potential Quadrature")
          .change_submod("Density on a grid", rho_mod);

        simde::type::tensor corr{{-1.83222, -0.402489},
                                 {-0.402489, -0.0886007}};
        auto error_mod =
          pluginplay::make_lambda<error_pt>([&](auto&&, auto&&, auto&& vxc_in) {
              REQUIRE(approximately_equal(vxc_in, corr, 1e-5));
              return uq_corr;
          });

        mod.change_submod("Incompleteness Error Model", error_mod);
        REQUIRE(mod.run_as<pt>(braket) == uq_corr);
    }
}
