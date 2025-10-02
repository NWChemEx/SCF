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

using grid_pt = simde::MolecularGrid;
using rho_pt  = simde::EDensityCollocationMatrix;
using wf_type = simde::type::rscf_wf;
using pt      = simde::eval_braket<wf_type, simde::type::XC_e_type, wf_type>;
using tensorwrapper::operations::approximately_equal;

TEST_CASE("LibXCEnergy") {
    using float_type = double;
    using index_set  = typename wf_type::orbital_index_set_type;

    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("LibXC Energy");

    auto rho    = test_scf::h2_density<float_type>();
    auto aos    = test_scf::h2_aos();
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

    auto exc = mod.run_as<pt>(braket);
    simde::type::tensor corr(-0.1031596741234082);

    REQUIRE(approximately_equal(exc, corr, 1e-6));
}
