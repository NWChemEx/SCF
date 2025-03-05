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

using rscf_wf = simde::type::rscf_wf;
using pt      = simde::InitialGuess<rscf_wf>;
using simde::type::tensor;

TEMPLATE_LIST_TEST_CASE("Core", "", test_scf::float_types) {
    using float_type = TestType;

    auto mm  = test_scf::load_modules<float_type>();
    auto aos = test_scf::h2_aos();
    auto H   = test_scf::h2_hamiltonian();

    auto mod = mm.at("Core guess");
    auto psi = mod.template run_as<pt>(H, aos);

    typename rscf_wf::orbital_index_set_type occs{0};
    REQUIRE(psi.orbital_indices() == occs);
    REQUIRE(psi.orbitals().from_space() == aos);
    const auto& evals       = psi.orbitals().diagonalized_matrix();
    using allocator_type    = tensorwrapper::allocator::Eigen<float_type>;
    const auto& eval_buffer = allocator_type::rebind(evals.buffer());

    const auto tol = 1E-6;
    using Catch::Matchers::WithinAbs;
    if constexpr(std::is_same_v<float_type, double>) {
        REQUIRE_THAT(eval_buffer.at(0), WithinAbs(-1.25330893, tol));
        REQUIRE_THAT(eval_buffer.at(1), WithinAbs(-0.47506974, tol));
    }
}