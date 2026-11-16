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

using rho_pt = simde::EDensityCollocationMatrix;
using pt     = scf::property_types::xc_energy_quadrature;
using tensorwrapper::operations::approximately_equal;

TEST_CASE("LibXCEnergyQuadrature") {
    using float_type = double;

    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("LibXC Energy Quadrature");

    auto rho    = test_scf::h2_density<float_type>();
    auto psi    = test_scf::h2_wave_function<float_type>();
    auto func   = chemist::qm_operator::xc_functional::SVWN3;
    auto xc_hat = test_scf::h2_xc<float_type>(func);
    chemist::braket::BraKet braket(psi, xc_hat, psi);

    // I compared these values from a modified Psi4numpy example with 6 angular
    // and 2 radial points per atom using Psi4's "SVWN" functional. These are
    // the points for the first batch. The resulting E_xc disagrees with Psi4 by
    // 1e-2; however, using the "correct" values from libxc_lda_energy_density
    // (which differ from Psi4 by about 1e-3) reproduces the correct value here
    //  to 1e-6, i.e., this test is consistent with the libxc_lda_energy_density
    // test.
    auto grid = test_scf::h2_xc_grid();

    auto rho_mod =
      pluginplay::make_lambda<rho_pt>([&](auto&& grid_in, auto&& P) {
          REQUIRE(grid_in == grid);
          REQUIRE(P == rho);
          return test_scf::h2_rho_on_xc_grid();
      });

    mod.change_submod("Density on a grid", rho_mod);

#ifdef BUILD_LIBXC
    auto exc = mod.run_as<pt>(braket, grid);
    simde::type::tensor corr(-0.1031596741234082);
    REQUIRE(approximately_equal(exc, corr, 1e-6));
#else
    REQUIRE_THROWS_AS(mod.run_as<pt>(braket, grid), std::runtime_error);
#endif
}
