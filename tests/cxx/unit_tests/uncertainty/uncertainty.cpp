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

#include "../../test_scf.hpp"
#include <vector>

using tensor_t   = simde::type::tensor;
using mo_pt      = scf::property_types::MOUncertainty;
using conv_pt    = scf::property_types::ConvergenceUncertainty;
using shape_type = tensorwrapper::shape::Smooth;

namespace {

/** Makes a tensor of FloatType from centers and radii.
 *
 *  The radii are ignored if FloatType is not an uncertain type.
 */
template<typename FloatType>
tensor_t make_tensor(shape_type shape, const std::vector<double>& centers,
                     const std::vector<double>& radii = {}) {
    using tensorwrapper::types::construct_uq_type;
    std::vector<FloatType> data;
    for(std::size_t i = 0; i < centers.size(); ++i) {
        const auto r = radii.empty() ? 0.0 : radii[i];
        if constexpr(tensorwrapper::types::is_uq_type_v<FloatType>)
            data.push_back(construct_uq_type<FloatType>(centers[i], r));
        else
            data.push_back(centers[i]);
    }
    tensorwrapper::buffer::Contiguous buffer(std::move(data), shape);
    return tensor_t(shape, std::move(buffer));
}

/// The centers and radii of the elements of a tensor of FloatType
template<typename FloatType>
auto centers_and_radii(const tensor_t& t) {
    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    using tensorwrapper::types::uq_center;
    using tensorwrapper::types::uq_upper;
    std::vector<double> centers, radii;
    for(const auto& x : get_raw_data<FloatType>(make_contiguous(t.buffer()))) {
        centers.push_back(double(uq_center(x)));
        radii.push_back(double(uq_upper(x) - uq_center(x)));
    }
    return std::pair{centers, radii};
}

// A two-orbital problem with S = 1. C and eps are the exact eigenpairs of
// the center of F
const shape_type matrix{2, 2};
const shape_type vector{2};
const shape_type scalar{};
const std::vector<double> F_centers{1.0, 0.0, 0.0, 2.0};
const std::vector<double> C_centers{1.0, 0.0, 0.0, 1.0};
const std::vector<double> eps_centers{1.0, 2.0};

} // namespace

TEST_CASE("No MO uncertainty") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("No MO uncertainty");

    auto F   = make_tensor<double>(matrix, F_centers);
    auto S   = make_tensor<double>(matrix, C_centers);
    auto eps = make_tensor<double>(vector, eps_centers);
    auto C   = make_tensor<double>(matrix, C_centers);
    auto dP  = make_tensor<double>(matrix, {1.0E-4, 0.0, 0.0, 1.0E-4});

    REQUIRE(mod.run_as<mo_pt>(F, S, eps, C, dP) == C);
}

TEST_CASE("First-order MO uncertainty") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("First-order MO uncertainty");

    SECTION("double") {
        auto F   = make_tensor<double>(matrix, F_centers);
        auto S   = make_tensor<double>(matrix, C_centers);
        auto eps = make_tensor<double>(vector, eps_centers);
        auto C   = make_tensor<double>(matrix, C_centers);
        auto dP  = make_tensor<double>(matrix, {1.0E-4, 0.0, 0.0, 1.0E-4});

        // Plain doubles carry no uncertainty, so C is returned unchanged
        REQUIRE(mod.run_as<mo_pt>(F, S, eps, C, dP) == C);
    }

#ifdef ENABLE_SIGMA
    using idouble = tensorwrapper::types::idouble;

    auto S   = make_tensor<idouble>(matrix, C_centers);
    auto eps = make_tensor<idouble>(vector, eps_centers);
    auto C   = make_tensor<idouble>(matrix, C_centers);

    SECTION("Uncertainty from the Fock matrix") {
        const double r = 1.0E-3;
        auto F  = make_tensor<idouble>(matrix, F_centers, {0.0, r, r, 0.0});
        auto dP = make_tensor<idouble>(matrix, {0.0, 0.0, 0.0, 0.0});

        auto rv                      = mod.run_as<mo_pt>(F, S, eps, C, dP);
        const auto& [centers, radii] = centers_and_radii<idouble>(rv);

        // The coupling through F, divided by the gap of 1, mixes the orbitals
        for(std::size_t i = 0; i < 4; ++i)
            REQUIRE_THAT(centers[i],
                         Catch::Matchers::WithinAbs(C_centers[i], 1.0E-12));
        REQUIRE(radii[0] < 1.0E-12);
        REQUIRE_THAT(radii[1], Catch::Matchers::WithinRel(r, 1.0E-6));
        REQUIRE_THAT(radii[2], Catch::Matchers::WithinRel(r, 1.0E-6));
        REQUIRE(radii[3] < 1.0E-12);
    }

    SECTION("Uncertainty from incomplete convergence") {
        auto F = make_tensor<idouble>(matrix, F_centers);
        auto dP =
          make_tensor<idouble>(matrix, {1.0E-4, -2.0E-4, -2.0E-4, 3.0E-4});

        auto rv                      = mod.run_as<mo_pt>(F, S, eps, C, dP);
        const auto& [centers, radii] = centers_and_radii<idouble>(rv);

        // Each coefficient is inflated by |dP| of the matching element
        const std::vector<double> corr{1.0E-4, 2.0E-4, 2.0E-4, 3.0E-4};
        for(std::size_t i = 0; i < 4; ++i) {
            REQUIRE_THAT(centers[i],
                         Catch::Matchers::WithinAbs(C_centers[i], 1.0E-12));
            REQUIRE_THAT(radii[i], Catch::Matchers::WithinRel(corr[i], 1.0E-6));
        }
    }
#endif
}

TEST_CASE("No convergence uncertainty") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("No convergence uncertainty");

    auto E  = make_tensor<double>(scalar, {-1.0});
    auto dE = make_tensor<double>(scalar, {1.0E-5});
    REQUIRE(mod.run_as<conv_pt>(E, dE) == E);
}

TEST_CASE("Energy change uncertainty") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("Energy change uncertainty");

    SECTION("double") {
        auto E  = make_tensor<double>(scalar, {-1.0});
        auto dE = make_tensor<double>(scalar, {1.0E-5});
        REQUIRE(mod.run_as<conv_pt>(E, dE) == E);
    }

#ifdef ENABLE_SIGMA
    SECTION("interval") {
        using idouble = tensorwrapper::types::idouble;
        auto E        = make_tensor<idouble>(scalar, {-1.0}, {1.0E-6});
        // Only the center of the change is used, not its radius
        auto dE = make_tensor<idouble>(scalar, {-1.0E-5}, {1.0E-3});

        auto rv                      = mod.run_as<conv_pt>(E, dE);
        const auto& [centers, radii] = centers_and_radii<idouble>(rv);
        REQUIRE_THAT(centers[0], Catch::Matchers::WithinAbs(-1.0, 1.0E-12));
        REQUIRE_THAT(radii[0], Catch::Matchers::WithinRel(1.1E-5, 1.0E-6));
    }
#endif
}
