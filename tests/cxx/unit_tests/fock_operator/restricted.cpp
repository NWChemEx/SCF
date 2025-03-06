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

TEMPLATE_LIST_TEST_CASE("Restricted", "", test_scf::float_types) {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    using float_type   = double;
    using density_type = simde::type::decomposable_e_density;
    using pt           = simde::FockOperator<density_type>;
    auto H             = test_scf::h2_hamiltonian();
    auto h2            = test_scf::make_h2<simde::type::nuclei>();
    auto rho           = test_scf::h2_density<float_type>();

    SECTION("Many-electron version") {
        auto& mod = mm.at("Restricted Fock Op");

        SECTION("No density") {
            density_type rho_empty;
            auto F = mod.run_as<pt>(H, rho_empty);
            simde::type::many_electrons es(2);

            simde::type::fock F_corr;
            using simde::type::T_e_type;
            using simde::type::V_en_type;
            F_corr.emplace_back(1.0, std::make_unique<T_e_type>(es));
            F_corr.emplace_back(1.0, std::make_unique<V_en_type>(es, h2));
            REQUIRE(F == F_corr);
        }

        SECTION("Density") {
            auto F      = mod.run_as<pt>(H, rho);
            auto F_corr = test_scf::h2_fock<simde::type::many_electrons>();
            REQUIRE(F == F_corr);
        }
    }
    SECTION("One-electron version") {
        auto& mod = mm.at("Restricted One-Electron Fock Op");

        SECTION(" No density") {
            density_type rho_empty;
            auto f = mod.run_as<pt>(H, rho_empty);
            simde::type::fock f_corr;
            simde::type::electron e;
            using simde::type::t_e_type;
            using simde::type::v_en_type;
            f_corr.emplace_back(1.0, std::make_unique<t_e_type>(e));
            f_corr.emplace_back(1.0, std::make_unique<v_en_type>(e, h2));
            REQUIRE(f == f_corr);
        }

        SECTION("Density") {
            auto f      = mod.run_as<pt>(H, rho);
            auto f_corr = test_scf::h2_fock<simde::type::electron>();
            REQUIRE(f == f_corr);
        }
    }
}