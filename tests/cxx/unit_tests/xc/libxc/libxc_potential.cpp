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
using pt     = scf::property_types::xc_potential_quadrature;

TEST_CASE("LibXCPotentialQuadrature") {
    using float_type = double;
    using tensorwrapper::operations::approximately_equal;

    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("LibXC Potential Quadrature");

    auto rho  = test_scf::h2_density<float_type>();
    auto aos  = test_scf::h2_aos();
    auto func = chemist::qm_operator::xc_functional::SVWN3;

    // For this test only the weights matter (we hardcode the density array)
    auto grid = test_scf::h2_xc_grid();

    auto rho_mod =
      pluginplay::make_lambda<rho_pt>([&](auto&& grid_in, auto&& P) {
          REQUIRE(grid_in == grid);
          REQUIRE(P == rho);
          return test_scf::h2_rho_on_xc_grid();
      });

    mod.change_submod("Density on a grid", rho_mod);

    simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
    chemist::braket::BraKet braket(aos, xc_op, aos);

#ifdef BUILD_LIBXC
    auto vxc = mod.run_as<pt>(braket, grid);
    simde::type::tensor corr{{-1.83222, -0.402489}, {-0.402489, -0.0886007}};
    REQUIRE(approximately_equal(vxc, corr, 1e-5));
#else
    REQUIRE_THROWS_AS(mod.run_as<pt>(braket, grid), std::runtime_error);
#endif
}
