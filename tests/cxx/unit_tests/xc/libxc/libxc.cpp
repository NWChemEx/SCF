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
#include "xc.h"
#include "xc_funcs.h"
#include <scf/xc/libxc/libxc.hpp>

TEST_CASE("to_libxc_codes") {
    using namespace chemist::qm_operator;
    using namespace scf::xc::libxc;
#ifdef BUILD_LIBXC
    auto [x, c] = to_libxc_codes(xc_functional::SVWN3);
    REQUIRE(x == XC_LDA_X);
    REQUIRE(c == XC_LDA_C_VWN_3);

    auto [x2, c2] = to_libxc_codes(xc_functional::SVWN5);
    REQUIRE(x2 == XC_LDA_X);
    REQUIRE(c2 == XC_LDA_C_VWN);

    REQUIRE_THROWS_AS(to_libxc_codes(xc_functional::B3LYP), std::runtime_error);
#else
    REQUIRE_THROWS_AS(to_libxc_codes(xc_functional::SVWN3), std::runtime_error);
#endif
}

TEST_CASE("libxc_lda_energy_density") {
    using namespace chemist::qm_operator;
    using namespace scf::xc::libxc;
    using tensorwrapper::operations::approximately_equal;

    // I compared these values from a modified Psi4numpy example with 6 angular
    // and 2 radial points per atom using Psi4's "SVWN" functional. These are
    // the points for the first batch. The agreement between these values and
    // the ones output by Psi4 is roughly 1e-3.

    simde::type::tensor rho_on_grid{0.01313525, 0.32141338, 0.01313525,
                                    0.3214133};
#ifdef BUILD_LIBXC
    auto exc = libxc_lda_energy_density(xc_functional::SVWN3, rho_on_grid);
    simde::type::tensor corr{-0.00280601, -0.182648, -0.00280601, -0.182648};
    REQUIRE(approximately_equal(exc, corr, 1e-6));
#else
    REQUIRE_THROWS_AS(
      libxc_lda_energy_density(xc_functional::SVWN3, rho_on_grid),
      std::runtime_error);
#endif
}

TEST_CASE("libxc_lda_energy_density_derivative") {
    using namespace chemist::qm_operator;
    using namespace scf::xc::libxc;
    using tensorwrapper::operations::approximately_equal;

    // I compared these values from a modified Psi4numpy example with 6 angular
    // and 2 radial points per atom using Psi4's "SVWN" functional. These are
    // the points for the first batch. The agreement between these values and
    // the ones output by Psi4 is roughly 1e-3.

    simde::type::tensor rho_on_grid{0.01313525, 0.32141338, 0.01313525,
                                    0.3214133};
#ifdef BUILD_LIBXC
    auto exc =
      libxc_lda_energy_density_derivative(xc_functional::SVWN3, rho_on_grid);
    simde::type::tensor corr{-0.278091, -0.744822, -0.278091, -0.744822};
    REQUIRE(approximately_equal(exc, corr, 1e-6));
#else
    REQUIRE_THROWS_AS(
      libxc_lda_energy_density_derivative(xc_functional::SVWN3, rho_on_grid),
      std::runtime_error);
#endif
}

TEST_CASE("tensorify_weights") {
    using namespace scf::xc::libxc;
    using tensorwrapper::operations::approximately_equal;

    std::vector<chemist::GridPoint> grid_points;
    for(std::size_t i = 0; i < 5; ++i)
        grid_points.push_back({static_cast<double>(i), 0.0, 0.0, 0.0});
    chemist::Grid grid(grid_points.begin(), grid_points.end());

    auto rv = parallelzone::runtime::RuntimeView();
    auto w  = tensorify_weights(grid, rv);
    simde::type::tensor corr{0.0, 1.0, 2.0, 3.0, 4.0};
    REQUIRE(approximately_equal(w, corr, 1e-8));
}

TEST_CASE("weight_a_matrix") {
    using namespace scf::xc::libxc;
    using tensorwrapper::operations::approximately_equal;

    simde::type::tensor w{1.0, 2.0, 3.0};
    simde::type::tensor b{
      {1.0, 4.0, 7.0}, {2.0, 5.0, 8.0}, {3.0, 6.0, 9.0}, {4.0, 10.0, 11.0}};
    auto rv = weight_a_matrix(w, b);
    simde::type::tensor corr{{1.0, 8.0, 21.0},
                             {2.0, 10.0, 24.0},
                             {3.0, 12.0, 27.0},
                             {4.0, 20.0, 33.0}};
    REQUIRE(approximately_equal(rv, corr, 1e-8));
}

TEST_CASE("batched_dot") {
    using namespace scf::xc::libxc;
    using tensorwrapper::operations::approximately_equal;

    simde::type::tensor b{
      {1.0, 4.0, 7.0}, {2.0, 5.0, 8.0}, {3.0, 6.0, 9.0}, {4.0, 10.0, 11.0}};

    SECTION("Sum row") {
        auto rv = batched_dot(b, b);
        simde::type::tensor corr{30.0, 177.0, 315.0};
        REQUIRE(approximately_equal(rv, corr, 1e-8));
    }
    SECTION("sum column") {
        auto rv = batched_dot(b, b, false);
        simde::type::tensor corr{66.0, 93.0, 126.0, 237.0};
        REQUIRE(approximately_equal(rv, corr, 1e-8));
    }
}
