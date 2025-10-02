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
#include <iostream>
#include <pluginplay/pluginplay.hpp>
#include <scf/xc/libxc/libxc.hpp>

using namespace scf;

using pt = simde::AOCollocationMatrix;

using tensorwrapper::operations::approximately_equal;

TEST_CASE("AOsOnGrid") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("AOs on a grid");

    SECTION("He STO-3G on a realistic grid") {
        // Assumes build directory is under root and we are running from the
        // build directory
        auto path = std::filesystem::current_path().parent_path();
        path += "/tests/he_grid.txt";
        mm.change_input("Grid From File", "Path to Grid File", path);
        using grid_pt = simde::MolecularGrid;
        auto he       = test_scf::make_he<chemist::Molecule>();
        auto grid     = mm.at("Grid From File").run_as<grid_pt>(he);
        auto ao_basis = test_scf::he_aos();
        auto rv       = mod.run_as<pt>(grid, ao_basis.ao_basis_set());
        auto runtime  = mm.get_runtime();
        auto weights  = scf::xc::libxc::tensorify_weights(grid, runtime);
        auto temp     = scf::xc::libxc::weight_a_matrix(weights, rv);
        auto norm     = scf::xc::libxc::batched_dot(temp, rv, false);
        typename simde::type::tensor::vector_il_type il{1.0};
        simde::type::tensor corr(il);
        REQUIRE(approximately_equal(norm, corr, 1e-6));
    }

    SECTION("H2 STO-3G on a realistic grid") {
        // Assumes build directory is under root and we are running from the
        // build directory
        auto path = std::filesystem::current_path().parent_path();
        path += "/tests/h2_grid.txt";
        mm.change_input("Grid From File", "Path to Grid File", path);
        using grid_pt = simde::MolecularGrid;
        auto h2       = test_scf::make_h2<chemist::Molecule>();
        auto grid     = mm.at("Grid From File").run_as<grid_pt>(h2);
        auto ao_basis = test_scf::h2_aos();
        auto rv       = mod.run_as<pt>(grid, ao_basis.ao_basis_set());
        auto runtime  = mm.get_runtime();
        auto weights  = scf::xc::libxc::tensorify_weights(grid, runtime);
        auto temp     = scf::xc::libxc::weight_a_matrix(weights, rv);
        auto norms    = scf::xc::libxc::batched_dot(temp, rv, false);
        typename simde::type::tensor::vector_il_type il{1.0, 1.0};
        simde::type::tensor corr(il);
        REQUIRE(approximately_equal(norms, corr, 1e-5));
    }
}
