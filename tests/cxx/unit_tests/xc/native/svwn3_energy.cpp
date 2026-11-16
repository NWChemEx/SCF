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

#include "../../../test_scf.hpp"

using rho_pt = simde::EDensityCollocationMatrix;
using pt     = scf::property_types::xc_energy_quadrature;
using tensorwrapper::operations::approximately_equal;

// Uses the same grid, density, and reference value as the LibXC test (see
// libxc/libxc_energy.cpp), since for a spin-unpolarized density the native
// SVWN3 is identical to LibXC's.
TEST_CASE("SVWN3EnergyQuadrature") {
    using float_type = double;

    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("SVWN3 Energy Quadrature");

    auto rho  = test_scf::h2_density<float_type>();
    auto psi  = test_scf::h2_wave_function<float_type>();
    auto grid = test_scf::h2_xc_grid();

    auto rho_mod =
      pluginplay::make_lambda<rho_pt>([&](auto&& grid_in, auto&& P) {
          REQUIRE(grid_in == grid);
          REQUIRE(P == rho);
          return test_scf::h2_rho_on_xc_grid();
      });

    mod.change_submod("Density on a grid", rho_mod);

    SECTION("SVWN3") {
        auto func   = chemist::qm_operator::xc_functional::SVWN3;
        auto xc_hat = test_scf::h2_xc<float_type>(func);
        chemist::braket::BraKet braket(psi, xc_hat, psi);

        auto exc = mod.run_as<pt>(braket, grid);
        simde::type::tensor corr(-0.1031596741234082);
        REQUIRE(approximately_equal(exc, corr, 1e-6));
    }

    SECTION("Throws for other functionals") {
        auto func   = chemist::qm_operator::xc_functional::SVWN5;
        auto xc_hat = test_scf::h2_xc<float_type>(func);
        chemist::braket::BraKet braket(psi, xc_hat, psi);
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket, grid), std::runtime_error);
    }
}
