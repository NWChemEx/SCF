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
#include <integrals/integrals.hpp>
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

    SECTION("He") {
        // Assumes build directory is under root and we are running from the
        // build directory
        auto path = std::filesystem::current_path().parent_path();
        path += "/tests/he_grid.txt";
        mm.change_input("Grid From File", "Path to Grid File", path);
        auto rho  = test_scf::he_density<float_type>();
        auto aos  = test_scf::he_aos();
        auto func = chemist::qm_operator::xc_functional::SVWN3;
        simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
        chemist::braket::BraKet braket(aos, xc_op, aos);

#ifdef BUILD_LIBXC
        auto vxc = mod.run_as<pt>(braket);
        typename simde::type::tensor::matrix_il_type il{{-0.668319}};
        simde::type::tensor corr(il);
        REQUIRE(approximately_equal(vxc, corr, 1e-5));
#else
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
#endif
    }

    SECTION("H2") {
        auto path = std::filesystem::current_path().parent_path();
        path += "/tests/h2_grid.txt";
        mm.change_input("Grid From File", "Path to Grid File", path);
        auto rho  = test_scf::h2_density<float_type>();
        auto aos  = test_scf::h2_aos();
        auto func = chemist::qm_operator::xc_functional::SVWN3;
        simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
        chemist::braket::BraKet braket(aos, xc_op, aos);

#ifdef BUILD_LIBXC
        auto vxc = mod.run_as<pt>(braket);
        simde::type::tensor corr{{-0.453301, -0.296985},
                                 {-0.296985, -0.453301}};
        REQUIRE(approximately_equal(vxc, corr, 1e-5));
#else
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
#endif
    }
}
