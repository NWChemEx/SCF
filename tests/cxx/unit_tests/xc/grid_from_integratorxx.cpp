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
#include <tensorwrapper/types/floating_point.hpp>

using namespace scf;

using pt = simde::MolecularGrid;

namespace {

// Small, fast-to-generate grid size shared by all sections below. 6 is a
// valid Lebedev-Laikov angular grid size (IntegratorXX only accepts exact
// supported sizes).
constexpr std::size_t radial_size  = 5;
constexpr std::size_t angular_size = 6;

pluginplay::ModuleManager make_mm() {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("Grid From IntegratorXX");
    mod.change_input("Radial Size", radial_size);
    mod.change_input("Angular Size", angular_size);
    return mm;
}

} // namespace

TEMPLATE_TEST_CASE("GridFromIntegratorXX", "[xc][molecular-grid]", float,
                   double) {
    using T = TestType;

    auto mm   = make_mm();
    auto& mod = mm.at("Grid From IntegratorXX");

    const auto float_type_str =
      std::is_same_v<T, float> ? std::string("float") : std::string("double");
    mod.change_input("Float Type", float_type_str);

    SECTION("He: point count and type") {
        auto he   = test_scf::make_he<chemist::Molecule>();
        auto grid = mod.run_as<pt>(he);

        REQUIRE(grid.size() == radial_size * angular_size);
        for(const auto& point : grid) {
            REQUIRE_NOTHROW(point.get_weight().template value<T>());
            REQUIRE_NOTHROW(point.get_x().template value<T>());
            REQUIRE_NOTHROW(point.get_y().template value<T>());
            REQUIRE_NOTHROW(point.get_z().template value<T>());
        }
    }

    SECTION("H2: point count scales with number of atoms") {
        auto h2   = test_scf::make_h2<chemist::Molecule>();
        auto grid = mod.run_as<pt>(h2);

        REQUIRE(grid.size() == 2 * radial_size * angular_size);
    }

    SECTION("Weights are all positive (single-atom, unpartitioned regions)") {
        auto he   = test_scf::make_he<chemist::Molecule>();
        auto grid = mod.run_as<pt>(he);

        for(const auto& point : grid) {
            REQUIRE(point.get_weight().template value<T>() > T(0));
        }
    }
}

TEST_CASE("GridFromIntegratorXX: Radial Quadrature Type") {
    // Exercises every accepted "Radial Quadrature Type" string (parsed via
    // IntegratorXX::radial_from_string, not a local re-implementation), to
    // check they all resolve to a working grid rather than just the default
    // ("MuraKnowles") exercised elsewhere.
    auto he = test_scf::make_he<chemist::Molecule>();

    for(const auto& name :
        {"Becke", "MurrayHandyLaming", "MuraKnowles", "TreutlerAhlrichs"}) {
        // run_as() locks the module, so each name needs its own instance.
        auto mm   = make_mm();
        auto& mod = mm.at("Grid From IntegratorXX");
        mod.change_input("Radial Quadrature Type", std::string(name));
        auto grid = mod.run_as<pt>(he);
        REQUIRE(grid.size() == radial_size * angular_size);
    }
}

TEST_CASE("GridFromIntegratorXX: Float Type parsing") {
    auto mm   = make_mm();
    auto& mod = mm.at("Grid From IntegratorXX");
    auto he   = test_scf::make_he<chemist::Molecule>();

    SECTION("Invalid Float Type string") {
        mod.change_input("Float Type", std::string("not-a-real-type"));
        REQUIRE_THROWS_AS(mod.run_as<pt>(he), std::runtime_error);
    }

    SECTION("Invalid Pruning Scheme string") {
        mod.change_input("Pruning Scheme", std::string("not-a-real-scheme"));
        REQUIRE_THROWS_AS(mod.run_as<pt>(he), std::runtime_error);
    }

    SECTION("Invalid Partition Scheme string") {
        mod.change_input("Partition Scheme", std::string("not-a-real-scheme"));
        REQUIRE_THROWS_AS(mod.run_as<pt>(he), std::runtime_error);
    }

    SECTION("udouble Float Type without Sigma support") {
#ifndef ENABLE_SIGMA
        mod.change_input("Float Type", std::string("udouble"));
        REQUIRE_THROWS_AS(mod.run_as<pt>(he), std::runtime_error);
#endif
    }

#ifndef ENABLE_SIGMA
    SECTION("idouble Float Type without Sigma support") {
        mod.change_input("Float Type", std::string("idouble"));
        REQUIRE_THROWS_AS(mod.run_as<pt>(he), std::runtime_error);
    }
#endif

    // adouble/tadouble/tmdouble are not yet supported at all -- see the
    // case-by-case comment in grid_from_integratorxx.cpp -- and throw
    // unconditionally, regardless of ENABLE_SIGMA.
    SECTION("adouble Float Type: not yet supported") {
        mod.change_input("Float Type", std::string("adouble"));
        REQUIRE_THROWS_AS(mod.run_as<pt>(he), std::runtime_error);
    }

    SECTION("tadouble Float Type: not yet supported") {
        mod.change_input("Float Type", std::string("tadouble"));
        REQUIRE_THROWS_AS(mod.run_as<pt>(he), std::runtime_error);
    }

    SECTION("tmdouble Float Type: not yet supported") {
        mod.change_input("Float Type", std::string("tmdouble"));
        REQUIRE_THROWS_AS(mod.run_as<pt>(he), std::runtime_error);
    }
}

#ifdef ENABLE_SIGMA
TEST_CASE("GridFromIntegratorXX: udouble") {
    using T = tensorwrapper::types::udouble;

    auto mm   = make_mm();
    auto& mod = mm.at("Grid From IntegratorXX");
    mod.change_input("Float Type", std::string("udouble"));

    auto he   = test_scf::make_he<chemist::Molecule>();
    auto grid = mod.run_as<pt>(he);

    REQUIRE(grid.size() == radial_size * angular_size);
    for(const auto& point : grid) {
        REQUIRE_NOTHROW(point.get_weight().template value<T>());
        REQUIRE_NOTHROW(point.get_x().template value<T>());
        REQUIRE_NOTHROW(point.get_y().template value<T>());
        REQUIRE_NOTHROW(point.get_z().template value<T>());
    }
}

TEST_CASE("GridFromIntegratorXX: idouble") {
    using T = tensorwrapper::types::idouble;

    auto mm   = make_mm();
    auto& mod = mm.at("Grid From IntegratorXX");
    mod.change_input("Float Type", std::string("idouble"));

    SECTION("He: point count and type") {
        auto he   = test_scf::make_he<chemist::Molecule>();
        auto grid = mod.run_as<pt>(he);

        REQUIRE(grid.size() == radial_size * angular_size);
        for(const auto& point : grid) {
            REQUIRE_NOTHROW(point.get_weight().template value<T>());
            REQUIRE_NOTHROW(point.get_x().template value<T>());
            REQUIRE_NOTHROW(point.get_y().template value<T>());
            REQUIRE_NOTHROW(point.get_z().template value<T>());
        }
    }

    SECTION("H2: point count scales with number of atoms, weights are "
            "partitioned (non-negative, single-atom regions)") {
        auto h2   = test_scf::make_h2<chemist::Molecule>();
        auto grid = mod.run_as<pt>(h2);

        REQUIRE(grid.size() == 2 * radial_size * angular_size);
        for(const auto& point : grid) {
            // Checks the enclosure's actual lower bound, not just sigma's
            // (median-based) relational operators, so this is a genuine
            // non-negativity check rather than an artifact of that ordering.
            // Not `> 0.0`: make_h2's nuclei sit on the z-axis, and this
            // angular_size (6) is a Lebedev-Laikov grid with points exactly
            // at (0,0,+-1) -- i.e. exactly on the bond axis. For a point
            // beyond the bond midpoint along that direction, the true Becke
            // weight for the near atom is exactly 0, not just small, so a
            // sound/tight enclosure legitimately touches zero (observed as
            // -0.0) rather than staying strictly positive.
            auto weight = point.get_weight().template value<T>();
            REQUIRE(tensorwrapper::types::uq_lower(weight) >= 0.0);
        }
    }
}
#endif
