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

using pt = simde::aos_op_base_aos;
using simde::type::tensor;

using erased_type =
  chemist::braket::BraKet<simde::type::aos, simde::type::op_base_type,
                          simde::type::aos>;

TEMPLATE_LIST_TEST_CASE("SCFIntegralsDriver", "", test_scf::float_types) {
    using float_type = TestType;
    auto mm          = test_scf::load_modules<float_type>();
    auto aos         = test_scf::h2_aos();
    auto mod         = mm.at("SCF integral driver");
    simde::type::electron e;
    auto rho = test_scf::h2_density<float_type>();

    using tensorwrapper::operations::approximately_equal;

    SECTION("Calling Kinetic") {
        auto& tmod = mm.at("Kinetic");
        simde::type::t_e_type t_e(e);
        chemist::braket::BraKet braket(aos, t_e, aos);
        erased_type copy_braket(braket);
        const auto& T      = mod.template run_as<pt>(copy_braket);
        const auto& T_corr = tmod.template run_as<simde::aos_t_e_aos>(braket);
        REQUIRE(approximately_equal(T, T_corr, 1E-6));
    }

    SECTION("Calling Electron-Nuclear Attraction") {
        auto& tmod     = mm.at("Nuclear");
        auto h2_nuclei = test_scf::make_h2<simde::type::nuclei>();
        simde::type::v_en_type v_en(e, h2_nuclei);
        chemist::braket::BraKet braket(aos, v_en, aos);
        erased_type copy_braket(braket);
        const auto& V      = mod.template run_as<pt>(copy_braket);
        const auto& V_corr = tmod.template run_as<simde::aos_v_en_aos>(braket);
        REQUIRE(approximately_equal(V, V_corr, 1E-6));
    }

    SECTION("Calling J Matrix") {
        auto& jmod = mm.at("Four center J builder");
        simde::type::j_e_type j_e(e, rho);
        chemist::braket::BraKet braket(aos, j_e, aos);
        erased_type copy_braket(braket);
        const auto& J      = mod.template run_as<pt>(copy_braket);
        const auto& J_corr = jmod.template run_as<simde::aos_j_e_aos>(braket);
        REQUIRE(approximately_equal(J, J_corr, 1E-6));
    }

    SECTION("Calling K Matrix") {
        auto& kmod = mm.at("Four center K builder");
        simde::type::k_e_type k_e(e, rho);
        chemist::braket::BraKet braket(aos, k_e, aos);
        erased_type copy_braket(braket);
        const auto& K      = mod.template run_as<pt>(copy_braket);
        const auto& K_corr = kmod.template run_as<simde::aos_k_e_aos>(braket);
        REQUIRE(approximately_equal(K, K_corr, 1E-6));
    }

    SECTION("Calling density matrix") {
        auto& pmod = mm.at("Density matrix builder");
        auto cmos  = test_scf::h2_cmos<float_type>();
        std::vector<int> occs{1, 0};
        simde::type::rho_e<simde::type::cmos> rho_hat(cmos, occs);
        chemist::braket::BraKet braket(aos, rho_hat, aos);
        erased_type copy_braket(braket);
        const auto& P      = mod.template run_as<pt>(copy_braket);
        using op_pt        = simde::aos_rho_e_aos<simde::type::cmos>;
        const auto& P_corr = pmod.template run_as<op_pt>(braket);
        REQUIRE(approximately_equal(P, P_corr, 1E-6));
    }

    // SECTION("Calling Fock Matrix") {
    //     auto& fmod = mm.at("Fock matrix builder");
    //     auto f_e   = test_scf::h2_fock<simde::type::electron>();
    //     chemist::braket::BraKet braket(aos, f_e, aos);
    //     erased_type copy_braket(braket);
    //     const auto& F      = mod.run_as<pt>(copy_braket);
    //     const auto& F_corr = fmod.run_as<simde::aos_f_e_aos>(braket);
    //     compare_matrices(F, F_corr);
    // }

    SECTION("Calling XC Potential (GauXC)") {
#ifdef BUILD_GAUXC
        auto func       = chemist::qm_operator::xc_functional::PBE0;
        auto rho        = test_scf::h2_density<double>();
        const auto& aos = rho.basis_set();
        simde::type::xc_e_type xc_op(func, e, rho);
        chemist::braket::BraKet xc_ij(aos, xc_op, aos);
        erased_type copy_braket(xc_ij);
        auto vxc = mod.template run_as<pt>(copy_braket);
        simde::type::tensor corr{{-0.357302, -0.23347}, {-0.23347, -0.357302}};
        REQUIRE(approximately_equal(vxc, corr, 1E-5));
#endif
    }

    SECTION("Calling XC Potential (LibXC)") {
#ifdef BUILD_LIBXC
        // Wired explicitly so this stays a LibXC test even when BUILD_GAUXC
        // is on a
        mm.change_submod("SCF integral driver", "XC Potential",
                         "LibXC Potential");
        // mod above is a copy taken before the rewiring, so re-fetch.
        auto& xc_mod    = mm.at("SCF integral driver");
        auto func       = chemist::qm_operator::xc_functional::SVWN3;
        auto rho        = test_scf::h2_density<double>();
        const auto& aos = rho.basis_set();
        simde::type::xc_e_type xc_op(func, e, rho);
        chemist::braket::BraKet xc_ij(aos, xc_op, aos);
        erased_type copy_braket(xc_ij);
        auto vxc = xc_mod.template run_as<pt>(copy_braket);

        // Psi4 1.11 VBase, STO-3G, {x: LDA_X, c: LDA_C_VWN_3}, 150x194
        // unpruned Becke grid, at this same density -- the value the
        // LibXCPotential test in tests/cxx/integration_tests/xc/ checks too.
        simde::type::tensor corr{{-0.453301, -0.296985},
                                 {-0.296985, -0.453301}};
        REQUIRE(approximately_equal(vxc, corr, 1E-5));
#endif
    }
}
