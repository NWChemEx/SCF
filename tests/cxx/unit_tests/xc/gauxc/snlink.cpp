/*
 * Copyright 2022 NWChemEx-Project
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

using k_pt = simde::aos_k_e_aos;
using tensorwrapper::operations::approximately_equal;

TEST_CASE("snLinK") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("snLinK");

    SECTION("h2") {
        auto rho        = test_scf::h2_density<double>();
        const auto& aos = rho.basis_set();
        simde::type::electron e;
        simde::type::k_e_type k(e, rho);
        simde::type::braket k_ij(aos, k, aos);
        const auto& K = mod.run_as<k_pt>(k_ij);

        simde::type::tensor corr_k(
          {{0.627264, 0.561828}, {0.561828, 0.627264}});
        REQUIRE(approximately_equal(K, corr_k, 1E-5));
    }
}