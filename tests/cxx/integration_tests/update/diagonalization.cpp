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

using wf_t              = simde::type::rscf_wf;
using pt                = simde::UpdateGuess<wf_t>;
using orbital_index_set = typename wf_t::orbital_index_set_type;

using simde::type::t_e_type;
using simde::type::v_en_type;

TEST_CASE("Diagaonalization") {
    auto mm   = test_scf::load_modules();
    auto& mod = mm.at("Diagonalization Fock update");

    auto aos = test_scf::h2_aos();
    auto h2  = test_scf::make_h2<simde::type::nuclei>();
    simde::type::electron e;
    orbital_index_set occs{0};

    simde::type::fock f_e;
    f_e.emplace_back(1.0, std::make_unique<t_e_type>(e));
    f_e.emplace_back(1.0, std::make_unique<v_en_type>(e, h2));

    SECTION("No old guess") {
        simde::type::tensor empty;
        simde::type::cmos cmos(empty, aos, empty);
        simde::type::rscf_wf core_guess(occs, cmos);
        const auto& psi = mod.run_as<pt>(f_e, core_guess);
        REQUIRE(psi.orbital_indices() == occs);
        REQUIRE(psi.orbitals().from_space() == aos);
        const auto& evals       = psi.orbitals().diagonalized_matrix();
        using allocator_type    = tensorwrapper::allocator::Eigen<double, 1>;
        const auto& eval_buffer = allocator_type::rebind(evals.buffer());

        const auto tol = 1E-6;
        using Catch::Matchers::WithinAbs;
        REQUIRE_THAT(eval_buffer.value()(0), WithinAbs(-1.25330893, tol));
        REQUIRE_THAT(eval_buffer.value()(1), WithinAbs(-0.47506974, tol));
    }
}