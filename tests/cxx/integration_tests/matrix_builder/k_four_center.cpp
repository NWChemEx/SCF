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

using pt = simde::aos_k_e_aos;

TEST_CASE("KFourCenter") {
    auto mm   = test_scf::load_modules();
    auto& mod = mm.at("Four center K builder");
    auto aos  = test_scf::h2_aos();

    simde::type::k_e_type k_e(simde::type::electron{}, test_scf::h2_density());
    const auto& K = mod.run_as<pt>(chemist::braket::BraKet(aos, k_e, aos));

    using alloc_type    = tensorwrapper::allocator::Eigen<double, 2>;
    const auto& K_eigen = alloc_type::rebind(K.buffer());
    using Catch::Matchers::WithinAbs;
    REQUIRE_THAT(K_eigen.value()(0, 0), WithinAbs(0.627264, 1E-6));
    REQUIRE_THAT(K_eigen.value()(0, 1), WithinAbs(0.561828, 1E-6));
    REQUIRE_THAT(K_eigen.value()(1, 0), WithinAbs(0.561828, 1E-6));
    REQUIRE_THAT(K_eigen.value()(1, 1), WithinAbs(0.627264, 1E-6));
}