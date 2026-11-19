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

using xc_functional = chemist::qm_operator::xc_functional;

// The dispatcher should forward its inputs, unchanged, to the submodule for
// the operator's functional, and return that submodule's result unchanged.
TEST_CASE("XCDispatcher<aos_xc_e_aos>") {
    using pt = simde::aos_xc_e_aos;

    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("XC Potential Dispatcher");

    auto aos         = test_scf::h2_aos();
    auto rho         = test_scf::h2_density<double>();
    auto make_braket = [&](xc_functional func) {
        simde::type::xc_e_type xc_op(func, simde::type::electron{}, rho);
        return chemist::braket::BraKet(aos, xc_op, aos);
    };

    SECTION("Dispatches SVWN3") {
        auto braket = make_braket(xc_functional::SVWN3);
        simde::type::tensor corr{{1.0, 2.0}, {3.0, 4.0}};
        auto svwn3_mod = pluginplay::make_lambda<pt>([&](auto&& braket_in) {
            REQUIRE(braket_in == braket);
            return corr;
        });
        mod.change_submod("SVWN3", svwn3_mod);
        REQUIRE(mod.run_as<pt>(braket) == corr);
    }

    SECTION("Throws for an unsupported functional") {
        auto braket = make_braket(xc_functional::NONE);
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
    }

#ifndef BUILD_GAUXC
    SECTION("Throws for a functional needing GauXC") {
        auto braket = make_braket(xc_functional::PBE);
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
    }
#endif
}

TEST_CASE("XCDispatcher<xc_energy_pt>") {
    using wf_type   = simde::type::rscf_wf;
    using index_set = typename wf_type::orbital_index_set_type;
    using pt = simde::eval_braket<wf_type, simde::type::XC_e_type, wf_type>;

    pluginplay::ModuleManager mm;
    mm.load_plugin<scf::SCF>();
    auto& mod = mm.at("XC Energy Dispatcher");

    wf_type psi(index_set{0}, test_scf::h2_cmos<double>());
    auto make_braket = [&](xc_functional func) {
        return chemist::braket::BraKet(psi, test_scf::h2_xc<double>(func), psi);
    };

    SECTION("Dispatches SVWN3") {
        auto braket = make_braket(xc_functional::SVWN3);
        simde::type::tensor corr(-1.23);
        auto svwn3_mod = pluginplay::make_lambda<pt>([&](auto&& braket_in) {
            REQUIRE(braket_in == braket);
            return corr;
        });
        mod.change_submod("SVWN3", svwn3_mod);
        REQUIRE(mod.run_as<pt>(braket) == corr);
    }

    SECTION("Throws for an unsupported functional") {
        auto braket = make_braket(xc_functional::NONE);
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
    }

#ifndef BUILD_GAUXC
    SECTION("Throws for a functional needing GauXC") {
        auto braket = make_braket(xc_functional::PBE);
        REQUIRE_THROWS_AS(mod.run_as<pt>(braket), std::runtime_error);
    }
#endif
}
