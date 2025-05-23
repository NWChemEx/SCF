/*
 * Copyright 2024 NWChemEx-Project
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

#include "../test_scf.hpp"
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <scf/scf.hpp>
#include <simde/simde.hpp>

using pt = simde::charge_charge_interaction;
using Catch::Matchers::WithinAbs;

TEST_CASE("CoulombsLaw") {
    using float_type = double;
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("Coulomb's Law");

    tensorwrapper::allocator::Eigen<float_type> alloc(mm.get_runtime());
    tensorwrapper::shape::Smooth shape_corr{};
    auto pcorr = alloc.allocate(tensorwrapper::layout::Physical(shape_corr));
    using tensorwrapper::operations::approximately_equal;

    using charge_type = typename simde::type::charges::value_type;
    charge_type q0(1.0, 2.0, 3.0, 4.0);
    charge_type q1(-1.0, 3.0, 4.0, 5.0);
    charge_type q2(1.5, 4.0, 5.0, 6.0);

    simde::type::charges empty;
    simde::type::charges qs{q0, q1, q2};

    SECTION("empty points") {
        auto e = mod.run_as<pt>(empty, empty);
        pcorr->set_elem({}, 0.0);
        tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
        REQUIRE(approximately_equal(corr, e, 1E-6));
    }

    SECTION("charges w/ itself") {
        auto e = mod.run_as<pt>(qs, qs);
        pcorr->set_elem({}, -1.0103629710818451);
        tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
        REQUIRE(approximately_equal(corr, e, 1E-6));
    }

    SECTION("charges w/ empty") {
        auto e = mod.run_as<pt>(qs, empty);
        pcorr->set_elem({}, 0.0);
        tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
        REQUIRE(approximately_equal(corr, e, 1E-6));
    }

    SECTION("charges w/ different charges") {
        simde::type::charges qs0{q0};
        simde::type::charges qs12{q1, q2};
        auto e = mod.run_as<pt>(qs0, qs12);
        pcorr->set_elem({}, -0.1443375672974065);
        tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));
        REQUIRE(approximately_equal(corr, e, 1E-6));
    }
}