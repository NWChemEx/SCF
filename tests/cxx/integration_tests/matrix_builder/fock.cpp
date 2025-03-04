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

#include "../integration_tests.hpp"

using pt = simde::aos_f_e_aos;

using simde::type::t_e_type;
using simde::type::v_en_type;

TEST_CASE("Fock Matrix Builder") {
    auto mm   = test_scf::load_modules();
    auto& mod = mm.at("Fock Matrix Builder");
    auto aos  = test_scf::h2_aos();

    SECTION("No J or K") {
        auto h2 = test_scf::make_h2<simde::type::nuclei>();
        simde::type::electron e;

        simde::type::fock f_e;
        f_e.emplace_back(1.0, std::make_unique<t_e_type>(e));
        f_e.emplace_back(1.0, std::make_unique<v_en_type>(e, h2));
        const auto& F = mod.run_as<pt>(chemist::braket::BraKet(aos, f_e, aos));

        using alloc_type    = tensorwrapper::allocator::Eigen<double>;
        const auto& F_eigen = alloc_type::rebind(F.buffer());
        using Catch::Matchers::WithinAbs;
        REQUIRE_THAT(F_eigen.at(0, 0), WithinAbs(-1.120958, 1E-6));
        REQUIRE_THAT(F_eigen.at(0, 1), WithinAbs(-0.959374, 1E-6));
        REQUIRE_THAT(F_eigen.at(1, 0), WithinAbs(-0.959374, 1E-6));
        REQUIRE_THAT(F_eigen.at(1, 1), WithinAbs(-1.120958, 1E-6));
    }

    SECTION("With J and K") {
        auto f_e = test_scf::h2_fock<simde::type::electron>();

        const auto& F = mod.run_as<pt>(chemist::braket::BraKet(aos, f_e, aos));

        using alloc_type    = tensorwrapper::allocator::Eigen<double>;
        const auto& F_eigen = alloc_type::rebind(F.buffer());

        using Catch::Matchers::WithinAbs;
        REQUIRE_THAT(F_eigen.at(0, 0), WithinAbs(-0.319459, 1E-6));
        REQUIRE_THAT(F_eigen.at(0, 1), WithinAbs(-0.571781, 1E-6));
        REQUIRE_THAT(F_eigen.at(1, 0), WithinAbs(-0.571781, 1E-6));
        REQUIRE_THAT(F_eigen.at(1, 1), WithinAbs(-0.319459, 1E-6));
    }
}