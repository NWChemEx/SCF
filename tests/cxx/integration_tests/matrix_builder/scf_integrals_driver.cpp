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

namespace {

void compare_matrices(const tensor& A, const tensor& A_corr) {
    using Catch::Matchers::WithinAbs;
    using alloc_type          = tensorwrapper::allocator::Eigen<double, 2>;
    const auto& A_buffer      = alloc_type::rebind(A.buffer());
    const auto& A_corr_buffer = alloc_type::rebind(A_corr.buffer());
    const auto& A_eigen       = A_buffer.value();
    const auto& A_corr_eigen  = A_corr_buffer.value();

    const auto tol = 1E-6;

    REQUIRE_THAT(A_eigen(0, 0), WithinAbs(A_corr_eigen(0, 0), 1E-6));
    REQUIRE_THAT(A_eigen(0, 1), WithinAbs(A_corr_eigen(0, 1), 1E-6));
    REQUIRE_THAT(A_eigen(1, 0), WithinAbs(A_corr_eigen(1, 0), 1E-6));
    REQUIRE_THAT(A_eigen(1, 1), WithinAbs(A_corr_eigen(1, 1), 1E-6));
}

} // namespace

using erased_type =
  chemist::braket::BraKet<simde::type::aos, simde::type::op_base_type,
                          simde::type::aos>;

TEST_CASE("SCFIntegralsDriver") {
    auto mm  = test_scf::load_modules();
    auto aos = test_scf::h2_aos();
    auto mod = mm.at("SCF integral driver");
    simde::type::electron e;
    auto rho = test_scf::h2_density();

    SECTION("Calling Kinetic") {
        auto& tmod = mm.at("Kinetic");
        simde::type::t_e_type t_e(e);
        chemist::braket::BraKet braket(aos, t_e, aos);
        erased_type copy_braket(braket);
        const auto& T      = mod.run_as<pt>(copy_braket);
        const auto& T_corr = tmod.run_as<simde::aos_t_e_aos>(braket);
        compare_matrices(T, T_corr);
    }

    SECTION("Calling Electron-Nuclear Attraction") {
        auto& tmod     = mm.at("Nuclear");
        auto h2_nuclei = test_scf::make_h2<simde::type::nuclei>();
        simde::type::v_en_type v_en(e, h2_nuclei);
        chemist::braket::BraKet braket(aos, v_en, aos);
        erased_type copy_braket(braket);
        const auto& V      = mod.run_as<pt>(copy_braket);
        const auto& V_corr = tmod.run_as<simde::aos_v_en_aos>(braket);
        compare_matrices(V, V_corr);
    }

    SECTION("Calling J Matrix") {
        auto& jmod = mm.at("Four center J builder");
        simde::type::j_e_type j_e(e, rho);
        chemist::braket::BraKet braket(aos, j_e, aos);
        erased_type copy_braket(braket);
        const auto& J      = mod.run_as<pt>(copy_braket);
        const auto& J_corr = jmod.run_as<simde::aos_j_e_aos>(braket);
        compare_matrices(J, J_corr);
    }

    SECTION("Calling K Matrix") {
        auto& kmod = mm.at("Four center K builder");
        simde::type::k_e_type k_e(e, rho);
        chemist::braket::BraKet braket(aos, k_e, aos);
        erased_type copy_braket(braket);
        const auto& K      = mod.run_as<pt>(copy_braket);
        const auto& K_corr = kmod.run_as<simde::aos_k_e_aos>(braket);
        compare_matrices(K, K_corr);
    }

    SECTION("Calling density matrix") {
        auto& pmod = mm.at("Density matrix builder");
        auto cmos  = test_scf::h2_cmos();
        std::vector<int> occs{1, 0};
        simde::type::rho_e<simde::type::cmos> rho_hat(cmos, occs);
        chemist::braket::BraKet braket(aos, rho_hat, aos);
        erased_type copy_braket(braket);
        const auto& P      = mod.run_as<pt>(copy_braket);
        using op_pt        = simde::aos_rho_e_aos<simde::type::cmos>;
        const auto& P_corr = pmod.run_as<op_pt>(braket);
        compare_matrices(P, P_corr);
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
}