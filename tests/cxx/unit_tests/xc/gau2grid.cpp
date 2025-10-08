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

TEST_CASE("Gau2Grid") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("Gau2Grid");

    using float_type     = double;
    using prim_type      = chemist::basis_set::Primitive<float_type>;
    using cg_type        = chemist::basis_set::ContractedGaussian<prim_type>;
    using shell_type     = chemist::basis_set::Shell<cg_type>;
    using atomic_bs_type = chemist::basis_set::AtomicBasisSet<shell_type>;
    using ao_basis_type  = chemist::basis_set::AOBasisSet<atomic_bs_type>;

    // TODO: better testing auto cart = chemist::ShellType::cartesian;
    auto pure = chemist::ShellType::pure;

    SECTION("Manual examples") {
        // Taken from https://gau2grid.readthedocs.io/en/latest/c_example.html

        std::vector<chemist::GridPoint> grid_points;
        for(std::size_t i = 0; i < 5; ++i)
            grid_points.push_back({1.0, 0.0, 0.0, static_cast<double>(i)});
        chemist::Grid grid(grid_points.begin(), grid_points.end());

        ao_basis_type ao_basis;
        atomic_bs_type abs("n/a", 1, {0.0, 0.0, 0.0});

        std::vector<float_type> coefs{1.0};
        std::vector<float_type> exps{1.0};
        cg_type cg(coefs.begin(), coefs.end(), exps.begin(), exps.end(), 0.0,
                   0.0, 0.0);

        SECTION("Single Basis Function example") {
            abs.add_shell(pure, 0, cg);
            ao_basis.add_center(std::move(abs));
            auto rv = mod.run_as<pt>(grid, ao_basis);
            simde::type::tensor corr{
              {1.0, 0.367879, 0.0183156, 0.00012341, 0.0}};
            REQUIRE(approximately_equal(rv, corr, 1e-6));
        }

        SECTION("Multiple Basis Function example") {
            abs.add_shell(pure, 0, cg);
            abs.add_shell(pure, 1, cg);
            abs.add_shell(pure, 2, cg);
            ao_basis.add_center(std::move(abs));
            auto rv = mod.run_as<pt>(grid, ao_basis);
            simde::type::tensor corr{
              {1, 0.367879, 0.0183156, 0.00012341, 1.12535e-07},
              {0, 0, 0, 0, 0},
              {0, 0.367879, 0.0366313, 0.000370229, 4.50141e-07},
              {0, 0, 0, 0, 0},
              {0, 0, 0, 0, 0},
              {0, 0, 0, 0, 0},
              {0, 0.367879, 0.0732626, 0.00111069, 1.80056e-06},
              {0, 0, 0, 0, 0},
              {0, 0, 0, 0, 0}};
            REQUIRE(approximately_equal(rv, corr, 1e-6));
        };
    }

    SECTION("A single s Gaussian on a realistic grid") {
        // Assumes build directory is under root and we are running from the
        // build directory
        auto path = std::filesystem::current_path().parent_path();
        path += "/tests/he_grid.txt";
        mm.change_input("Grid From File", "Path to Grid File", path);
        using grid_pt = simde::MolecularGrid;
        auto he       = test_scf::make_he<chemist::Molecule>();
        auto grid     = mm.at("Grid From File").run_as<grid_pt>(he);
        ao_basis_type ao_basis;
        atomic_bs_type abs("n/a", 1, {0.0, 0.0, 0.0});

        std::vector<float_type> coefs{1.0};
        std::vector<float_type> exps{1.0};
        cg_type cg(coefs.begin(), coefs.end(), exps.begin(), exps.end(), 0.0,
                   0.0, 0.0);
        abs.add_shell(pure, 0, cg);
        ao_basis.add_center(std::move(abs));
        auto rv      = mod.run_as<pt>(grid, ao_basis);
        auto runtime = mm.get_runtime();
        auto weights = scf::xc::libxc::tensorify_weights(grid, runtime);
        simde::type::tensor sum;
        sum("m") = weights("i") * rv("m,i");
        // \int d^3r exp(-r^2) = 4pi \int_0^inf dr r^2 exp(-r^2) = pi^(3/2)
        simde::type::tensor corr{M_PI * std::sqrt(M_PI)};
        REQUIRE(approximately_equal(sum, corr, 1e-6));
    }
}
