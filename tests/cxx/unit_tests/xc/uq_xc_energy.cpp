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
using wf_type  = simde::type::rscf_wf;
using pt       = simde::eval_braket<wf_type, simde::type::XC_e_type, wf_type>;
using quad_pt  = scf::property_types::xc_energy_quadrature;
using error_pt = scf::property_types::xc_energy_uncertainty;
using tensorwrapper::operations::approximately_equal;

TEST_CASE("UQXCEnergy") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("UQ XC Energy");

    auto psi    = test_scf::h2_wave_function<double>();
    auto func   = chemist::qm_operator::xc_functional::SVWN3;
    auto xc_hat = test_scf::h2_xc<double>(func);
    chemist::braket::BraKet braket(psi, xc_hat, psi);
    auto grid = test_scf::h2_xc_grid();

    auto grid_mod = pluginplay::make_lambda<grid_pt>([&](auto&& mol) {
        REQUIRE(mol.size() == 2);
        return grid;
    });
    mod.change_submod("Integration grid", grid_mod);

    // Stands in for the uncertain result the error model returns
    simde::type::tensor uq_corr(4.56);

    SECTION("Forwards the grid and the quadrature to the error model") {
        simde::type::tensor exc(1.23);
        auto quad_mod = pluginplay::make_lambda<quad_pt>(
          [&](auto&& braket_in, auto&& grid_in) {
              REQUIRE(braket_in == braket);
              REQUIRE(grid_in == grid);
              return exc;
          });

        auto error_mod = pluginplay::make_lambda<error_pt>(
          [&](auto&& braket_in, auto&& grid_in, auto&& exc_in) {
              REQUIRE(braket_in == braket);
              REQUIRE(grid_in == grid);
              REQUIRE(exc_in == exc);
              return uq_corr;
          });

        mod.change_submod("XC Energy", quad_mod);
        mod.change_submod("Incompleteness Error Model", error_mod);
        REQUIRE(mod.run_as<pt>(braket) == uq_corr);
    }

    SECTION("Defaults to the native SVWN3 quadrature") {
        // Hardcode the density on the grid (see the SVWN3EnergyQuadrature
        // test, whose reference value this is)
        auto rho_mod = pluginplay::make_lambda<rho_pt>(
          [&](auto&&, auto&&) { return test_scf::h2_rho_on_xc_grid(); });
        mm.at("SVWN3 Energy Quadrature")
          .change_submod("Density on a grid", rho_mod);

        simde::type::tensor corr(-0.1031596741234082);
        auto error_mod =
          pluginplay::make_lambda<error_pt>([&](auto&&, auto&&, auto&& exc_in) {
              REQUIRE(approximately_equal(exc_in, corr, 1e-6));
              return uq_corr;
          });

        mod.change_submod("Incompleteness Error Model", error_mod);
        REQUIRE(mod.run_as<pt>(braket) == uq_corr);
    }
}
