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
#include <chemist/experimental/basis_set/contracted_gaussian.hpp>
#include <iostream>
#include <pluginplay/pluginplay.hpp>
#include <scf/xc/libxc/libxc.hpp>

using namespace scf;

using pt = simde::AOCollocationMatrix;

using tensorwrapper::operations::approximately_equal;

namespace {

using float_type     = double;
using prim_type      = chemist::basis_set::Primitive<float_type>;
using cg_type        = chemist::basis_set::ContractedGaussian<prim_type>;
using shell_type     = chemist::basis_set::Shell<cg_type>;
using atomic_bs_type = chemist::basis_set::AtomicBasisSet<shell_type>;
using ao_basis_type  = chemist::basis_set::AOBasisSet<atomic_bs_type>;

// Folds N^chi and N^G into the coefficients, as gau2grid expects
std::vector<double> fold_normalization(const std::vector<double>& cs,
                                       const std::vector<double>& zetas,
                                       unsigned l) {
    chemist::experimental::ContractedGaussian cg(
      cs.begin(), cs.end(), zetas.begin(), zetas.end(), l, 0.0, 0.0, 0.0);
    auto n_g = cg.normalization_constant();
    std::vector<double> rv(cs);
    for(std::size_t p = 0; p < cs.size(); ++p) {
        auto n_chi = cg[p].normalization_constant();
        rv[p] *=
          wtf::fp::float_cast<double>(n_g) * wtf::fp::float_cast<double>(n_chi);
    }
    return rv;
}

// Basis set with a single shell, centered on (x, y, z)
ao_basis_type one_shell_basis(chemist::ShellType purity, unsigned l,
                              const std::vector<double>& cs,
                              const std::vector<double>& zetas, double x,
                              double y, double z) {
    atomic_bs_type abs("n/a", 1, {x, y, z});
    cg_type cg(cs.begin(), cs.end(), zetas.begin(), zetas.end(), x, y, z);
    abs.add_shell(purity, l, cg);
    ao_basis_type ao_basis;
    ao_basis.add_center(std::move(abs));
    return ao_basis;
}

} // namespace

TEST_CASE("AOsOnGrid") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("AOs on a grid");

    auto path = test_scf::get_test_directory_path();

    SECTION("He STO-3G on a realistic grid") {
        path += "/he_grid.txt";
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
        path += "/h2_grid.txt";
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

    SECTION("Agrees with gau2grid") {
        // A single point, displaced from the center in all three directions
        std::vector<chemist::GridPoint> grid_points{{1.0, 0.3, -0.4, 0.5}};
        chemist::Grid grid(grid_points.begin(), grid_points.end());
        const double x = 0.1, y = 0.2, z = -0.3;

        // Several primitives so N^G is not 1
        std::vector<double> cs{0.15, 0.55, 0.45};
        std::vector<double> zetas{3.0, 1.0, 0.25};

        auto& g2g = mm.at("Gau2Grid");

        const auto pure = chemist::ShellType::pure;
        const auto cart = chemist::ShellType::cartesian;

        // Gau2Grid needs the normalization folded into the coefficients
        auto run_g2g = [&](chemist::ShellType purity, unsigned l) {
            auto norm_cs = fold_normalization(cs, zetas, l);
            auto bs      = one_shell_basis(purity, l, norm_cs, zetas, x, y, z);
            return g2g.run_as<pt>(grid, bs);
        };

        auto run_mod = [&](chemist::ShellType purity, unsigned l) {
            auto bs = one_shell_basis(purity, l, cs, zetas, x, y, z);
            return mod.run_as<pt>(grid, bs);
        };

        for(unsigned l = 0; l < 4; ++l) {
            DYNAMIC_SECTION("Spherical, l = " << l) {
                auto rv   = run_mod(pure, l);
                auto corr = run_g2g(pure, l);
                REQUIRE(approximately_equal(rv, corr, 1e-10));
            }
            DYNAMIC_SECTION("Cartesian, l = " << l) {
                auto rv   = run_mod(cart, l);
                auto corr = run_g2g(cart, l);
                REQUIRE(approximately_equal(rv, corr, 1e-10));
            }
        }

        SECTION("s is the same for both purities") {
            REQUIRE(
              approximately_equal(run_mod(pure, 0), run_mod(cart, 0), 1e-12));
        }

        SECTION("p is the same for both purities, up to ordering") {
            // CCA orders spherical p as m = -1, 0, 1, i.e., y, z, x
            auto sph    = run_mod(pure, 1);
            auto cart_p = run_mod(cart, 1);
            using tensorwrapper::buffer::make_contiguous;
            using wtf::fp::float_cast;
            const auto& sph_buf  = make_contiguous(sph.buffer());
            const auto& cart_buf = make_contiguous(cart_p.buffer());
            const std::vector<std::size_t> cart_row{1, 2, 0};
            for(std::size_t m = 0; m < 3; ++m) {
                auto lhs = float_cast<double>(sph_buf.get_elem({m, 0}));
                auto rhs =
                  float_cast<double>(cart_buf.get_elem({cart_row[m], 0}));
                REQUIRE(lhs == Catch::Approx(rhs).epsilon(1e-12));
            }
        }
    }
}
