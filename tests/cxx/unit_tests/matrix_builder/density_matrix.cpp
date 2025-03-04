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

#include "../../test_scf.hpp"

using pt = simde::aos_rho_e_aos<simde::type::cmos>;

TEST_CASE("Density Matrix Builder") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("Density matrix builder");
    auto aos  = test_scf::h2_aos();
    auto cmos = test_scf::h2_cmos();
    std::vector<int> occs{1, 0};
    simde::type::rho_e<simde::type::cmos> rho_hat(cmos, occs);

    chemist::braket::BraKet p_mn(aos, rho_hat, aos);
    const auto& P        = mod.run_as<pt>(p_mn);
    using allocator_type = tensorwrapper::allocator::Eigen<double>;
    const auto& P_eigen  = allocator_type::rebind(P.buffer());

    using Catch::Matchers::WithinAbs;
    REQUIRE_THAT(P_eigen.at(0, 0), WithinAbs(0.31980835, 1E-6));
    REQUIRE_THAT(P_eigen.at(0, 1), WithinAbs(0.31980835, 1E-6));
    REQUIRE_THAT(P_eigen.at(1, 0), WithinAbs(0.31980835, 1E-6));
    REQUIRE_THAT(P_eigen.at(1, 1), WithinAbs(0.31980835, 1E-6));
}