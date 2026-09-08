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
using rho_pt  = simde::EDensityCollocationMatrix;
using wf_type = simde::type::rscf_wf;
using pt      = simde::eval_braket<wf_type, simde::type::XC_e_type, wf_type>;

TEST_CASE("LibXCEnergy") {
    using float_type = double;
    using tensorwrapper::operations::approximately_equal;

    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    mm.change_submod("LibXC Energy", "Integration grid",
                     "Grid From IntegratorXX");
    auto& mod = mm.at("LibXC Energy");

    // Matches the (radial x angular) point count of the reference
    // he_grid.txt/h2_grid.txt files (150*194 = 29100, 2*150*194 = 58200)
    // generated with IntegratorXX's own defaults (MuraKnowles radial,
    // Lebedev-Laikov angular, unpruned, Becke partition), so this
    // cross-validates the generated grid against the known-good file-based
    // one -- not a bit-identical grid, hence the looser tolerance below
    // (vs. the GridFromFile-based unit tests' 1e-5).
    mm.change_input("Grid From IntegratorXX", "Radial Size", std::size_t(150));
    mm.change_input("Grid From IntegratorXX", "Angular Size", std::size_t(194));

    SECTION("He") {
        auto rho    = test_scf::he_density<float_type>();
        auto aos    = test_scf::he_aos();
        auto psi    = test_scf::he_wave_function<float_type>();
        auto func   = chemist::qm_operator::xc_functional::SVWN3;
        auto xc_hat = test_scf::he_xc<float_type>(func);
        chemist::braket::BraKet braket(psi, xc_hat, psi);

#ifdef BUILD_LIBXC
        auto exc = mod.run_as<pt>(braket);
        simde::type::tensor corr(-1.01982);
        REQUIRE(approximately_equal(exc, corr, 1e-3));
#else
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
#endif
    }

    SECTION("H2") {
        auto rho    = test_scf::h2_density<float_type>();
        auto aos    = test_scf::h2_aos();
        auto psi    = test_scf::h2_wave_function<float_type>();
        auto func   = chemist::qm_operator::xc_functional::SVWN3;
        auto xc_hat = test_scf::h2_xc<float_type>(func);
        chemist::braket::BraKet braket(psi, xc_hat, psi);

#ifdef BUILD_LIBXC
        auto exc = mod.run_as<pt>(braket);
        simde::type::tensor corr(-0.734458);
        REQUIRE(approximately_equal(exc, corr, 1e-3));
#else
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
#endif
    }
}
