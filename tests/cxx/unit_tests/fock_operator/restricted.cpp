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
#include <scf/scf.hpp>
#include <simde/simde.hpp>

TEST_CASE("Restricted") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);

    using density_type = simde::type::decomposable_e_density;
    using pt           = simde::FockOperator<density_type>;

    auto& mod = mm.at("Restricted Fock Op");

    SECTION("H2 Molecule") {
        auto H   = test_scf::h2_hamiltonian();
        auto rho = test_scf::h2_density();

        auto F      = mod.run_as<pt>(H, rho);
        auto F_corr = test_scf::h2_fock();
        REQUIRE(F == F_corr);
    }
}