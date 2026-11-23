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
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

using grid_pt   = simde::MolecularGrid;
using uq_pt     = simde::UQInitializer;
using e_quad_pt = scf::property_types::xc_energy_quadrature;
using v_quad_pt = scf::property_types::xc_potential_quadrature;
using simde::UQKind;
using simde::type::tensor;

namespace {

// Stands in for "Grid From IntegratorXX". It returns a one-point grid whose
// weight and x-coordinate are the radial and angular sizes it was run with,
// so the reference XC module can check them.
struct FakeGrid : pluginplay::ModuleBase {
    FakeGrid() : pluginplay::ModuleBase(this) {
        satisfies_property_type<grid_pt>();
        add_input<std::size_t>("Radial Size").set_default(std::size_t{0});
        add_input<std::size_t>("Angular Size").set_default(std::size_t{0});
    }

private:
    result_map_type run_(input_map_type inputs,
                         submodule_map_type) const override {
        const auto r = inputs.at("Radial Size").value<std::size_t>();
        const auto a = inputs.at("Angular Size").value<std::size_t>();
        auto rv      = results();
        return grid_pt::wrap_results(rv, sized_grid(r, a));
    }

public:
    /// The grid returned for a radial size @p r and angular size @p a
    static chemist::Grid sized_grid(std::size_t r, std::size_t a) {
        std::vector<chemist::GridPoint> points{chemist::GridPoint(
          static_cast<double>(r), static_cast<double>(a), 0.0, 0.0)};
        return chemist::Grid(points.begin(), points.end());
    }
};

/// The grid the given value was evaluated on, recording @p r and @p a
chemist::Grid given_grid(std::size_t r, std::size_t a) {
    auto grid = test_scf::h2_xc_grid();
    grid.set_radial_size(r);
    grid.set_angular_size(a);
    return grid;
}

/** @brief Runs a fresh copy of the error model @p key.
 *
 *  The given grid records the sizes @p r and @p a, and the reference grid
 *  must be built with the sizes @p ref_r and @p ref_a. The "Reference grid"
 *  submodule is a FakeGrid and the other submodules are lambdas.
 */
template<typename QuadPT, typename BraKetType>
tensor run_model(pluginplay::ModuleManager& mm, const std::string& key,
                 const BraKetType& braket, const tensor& value,
                 const tensor& ref, UQKind kind, const chemist::Grid& grid,
                 std::size_t ref_r, std::size_t ref_a) {
    using error_pt = scf::property_types::XCUncertain<QuadPT>;

    // Modules lock once they run, so each call uses its own copies
    static std::size_t counter = 0;
    const auto suffix          = " " + std::to_string(counter++);
    const auto copy_key        = key + suffix;
    const auto grid_key        = "Fake grid" + suffix;
    mm.copy_module(key, copy_key);
    mm.add_module<FakeGrid>(grid_key);
    auto& mod = mm.at(copy_key);

    const auto ref_grid = FakeGrid::sized_grid(ref_r, ref_a);
    mm.change_submod(copy_key, "Reference grid", grid_key);
    mod.change_submod(
      "Reference XC Module",
      pluginplay::make_lambda<QuadPT>([&](auto&& braket_in, auto&& grid_in) {
          REQUIRE(braket_in == braket);
          REQUIRE(grid_in == ref_grid);
          return ref;
      }));
    mod.change_submod("UQ Initializer", pluginplay::make_lambda<uq_pt>([=]() {
                          return simde::UQFactory(kind);
                      }));
    return mod.template run_as<error_pt>(braket, grid, value);
}

/// run_model for a given grid of 10 x 14, whose reference grid is 15 x 38
template<typename QuadPT, typename BraKetType>
tensor run_model(pluginplay::ModuleManager& mm, const std::string& key,
                 const BraKetType& braket, const tensor& value,
                 const tensor& ref, UQKind kind) {
    return run_model<QuadPT>(mm, key, braket, value, ref, kind,
                             given_grid(10, 14), 15, 38);
}

/// The elements of @p t (a scalar or a matrix) as T, in row-major order
template<typename T>
std::vector<T> elements(const tensor& t) {
    using tensorwrapper::buffer::make_contiguous;
    const auto& buffer = make_contiguous(t.buffer());
    const auto shape   = buffer.shape();
    std::vector<T> rv;
    if(shape.rank() == 0) {
        rv.push_back(wtf::fp::float_cast<T>(buffer.get_elem({})));
    } else {
        for(std::size_t i = 0; i < shape.extent(0); ++i)
            for(std::size_t j = 0; j < shape.extent(1); ++j)
                rv.push_back(wtf::fp::float_cast<T>(buffer.get_elem({i, j})));
    }
    return rv;
}

/// The elements of @p t, which holds doubles
std::vector<double> doubles(const tensor& t) { return elements<double>(t); }

#ifdef ENABLE_SIGMA
/// The (lower, nominal, upper) of each element of @p t, which holds UQType
template<typename UQType>
std::vector<std::array<double, 3>> bounds(const tensor& t) {
    using namespace tensorwrapper::types;
    std::vector<std::array<double, 3>> rv;
    for(const auto& x : elements<UQType>(t))
        rv.push_back({uq_lower(x), uq_center(x), uq_upper(x)});
    return rv;
}
#endif

} // namespace

TEST_CASE("ReferenceGridError (potential)") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    const auto key = "Reference Grid XC Potential Error";

    auto rho  = test_scf::h2_density<double>();
    auto aos  = test_scf::h2_aos();
    auto func = chemist::qm_operator::xc_functional::SVWN3;
    simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
    chemist::braket::BraKet braket(aos, xc_op, aos);

    // d = ref - value = {0.5, 0, -1, 0.25}: above, equal, below, above
    tensor value{{1.0, 2.0}, {3.0, 4.0}};
    tensor ref{{1.5, 2.0}, {2.0, 4.25}};

    SECTION("none returns the given value") {
        auto rv =
          run_model<v_quad_pt>(mm, key, braket, value, ref, UQKind::none);
        REQUIRE(doubles(rv) == doubles(value));
    }

    SECTION("throws if the value does not hold doubles") {
        tensorwrapper::shape::Smooth shape{2, 2};
        tensorwrapper::buffer::Contiguous buffer(
          std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f}, shape);
        tensor float_value(shape, std::move(buffer));
        REQUIRE_THROWS_AS(
          run_model<v_quad_pt>(mm, key, braket, float_value, ref, UQKind::none),
          std::runtime_error);
    }

    SECTION("reference grid sizes") {
        auto run = [&](std::size_t r, std::size_t a, std::size_t ref_r,
                       std::size_t ref_a) {
            return run_model<v_quad_pt>(mm, key, braket, value, ref,
                                        UQKind::none, given_grid(r, a), ref_r,
                                        ref_a);
        };
        // Radial: x 1.5, rounded up. Angular: the smallest Lebedev-Laikov
        // grid with at least twice as many points.
        REQUIRE(doubles(run(150, 194, 225, 434)) == doubles(value));
        REQUIRE(doubles(run(75, 302, 113, 770)) == doubles(value));
        // Angular: the largest grid if no grid has twice as many points
        REQUIRE(doubles(run(150, 3074, 225, 5810)) == doubles(value));
        // Angular: the largest grid is not refined
        REQUIRE(doubles(run(150, 5810, 225, 5810)) == doubles(value));
    }

    SECTION("throws if the given grid does not record its sizes") {
        auto no_sizes = test_scf::h2_xc_grid();
        REQUIRE_THROWS_AS(run_model<v_quad_pt>(mm, key, braket, value, ref,
                                               UQKind::none, no_sizes, 15, 38),
                          std::invalid_argument);

        auto radial_only = test_scf::h2_xc_grid();
        radial_only.set_radial_size(10);
        REQUIRE_THROWS_AS(run_model<v_quad_pt>(mm, key, braket, value, ref,
                                               UQKind::none, radial_only, 15,
                                               38),
                          std::invalid_argument);
    }

    SECTION("throws if the angular size is not a Lebedev-Laikov size") {
        REQUIRE_THROWS_AS(run_model<v_quad_pt>(mm, key, braket, value, ref,
                                               UQKind::none, given_grid(10, 15),
                                               15, 38),
                          std::invalid_argument);
    }

    SECTION("throws if the reference grid module has no size inputs") {
        const auto copy_key = std::string(key) + " no size inputs";
        mm.copy_module(key, copy_key);
        auto& mod = mm.at(copy_key);
        mod.change_submod("Reference grid",
                          pluginplay::make_lambda<grid_pt>(
                            [](auto&&) { return test_scf::h2_xc_grid(); }));
        mod.change_submod("Reference XC Module",
                          pluginplay::make_lambda<v_quad_pt>(
                            [&](auto&&, auto&&) { return ref; }));
        mod.change_submod("UQ Initializer",
                          pluginplay::make_lambda<uq_pt>(
                            []() { return simde::UQFactory(UQKind::none); }));
        using error_pt = scf::property_types::XCUncertain<v_quad_pt>;
        REQUIRE_THROWS_AS(
          mod.run_as<error_pt>(braket, given_grid(10, 14), value),
          std::invalid_argument);
    }

#ifdef ENABLE_SIGMA
    using tensorwrapper::types::idouble;
    using tensorwrapper::types::tpdouble;
    using triple = std::array<double, 3>;

    SECTION("triplex is symmetric") {
        auto rv =
          run_model<v_quad_pt>(mm, key, braket, value, ref, UQKind::triplex);
        std::vector<triple> corr{
          {0.5, 1.0, 1.5}, {2.0, 2.0, 2.0}, {2.0, 3.0, 4.0}, {3.75, 4.0, 4.25}};
        REQUIRE(bounds<tpdouble>(rv) == corr);
    }

    SECTION("interval keeps the bounds") {
        auto rv =
          run_model<v_quad_pt>(mm, key, braket, value, ref, UQKind::interval);
        auto b = bounds<idouble>(rv);
        std::vector<std::array<double, 2>> corr{
          {0.5, 1.5}, {2.0, 2.0}, {2.0, 4.0}, {3.75, 4.25}};
        for(std::size_t i = 0; i < corr.size(); ++i) {
            REQUIRE(b[i][0] == corr[i][0]);
            REQUIRE(b[i][2] == corr[i][1]);
        }
    }
#endif
}

TEST_CASE("ReferenceGridError (energy)") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    const auto key = "Reference Grid XC Energy Error";

    auto psi    = test_scf::h2_wave_function<double>();
    auto func   = chemist::qm_operator::xc_functional::SVWN3;
    auto xc_hat = test_scf::h2_xc<double>(func);
    chemist::braket::BraKet braket(psi, xc_hat, psi);
    tensor value(-1.0);

    SECTION("none returns the given value") {
        auto rv = run_model<e_quad_pt>(mm, key, braket, value, tensor(-0.5),
                                       UQKind::none);
        REQUIRE(doubles(rv) == doubles(value));
    }

#ifdef ENABLE_SIGMA
    using tensorwrapper::types::tpdouble;
    using triple = std::array<double, 3>;

    SECTION("triplex, reference above") {
        auto rv = run_model<e_quad_pt>(mm, key, braket, value, tensor(-0.5),
                                       UQKind::triplex);
        REQUIRE(bounds<tpdouble>(rv) ==
                std::vector<triple>{{-1.5, -1.0, -0.5}});
    }

    SECTION("triplex, reference below") {
        auto rv = run_model<e_quad_pt>(mm, key, braket, value, tensor(-1.5),
                                       UQKind::triplex);
        REQUIRE(bounds<tpdouble>(rv) ==
                std::vector<triple>{{-1.5, -1.0, -0.5}});
    }
#endif
}
