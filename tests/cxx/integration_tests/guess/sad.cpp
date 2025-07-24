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

using simde::type::tensor;
using shape_type   = tensorwrapper::shape::Smooth;
using cmos_type    = simde::type::cmos;
using density_type = simde::type::decomposable_e_density;
using rscf_wf      = simde::type::rscf_wf;
using occ_index    = typename rscf_wf::orbital_index_set_type;

using pt             = simde::InitialGuess<rscf_wf>;
using initial_rho_pt = simde::InitialDensity;

using tensorwrapper::operations::approximately_equal;

TEMPLATE_LIST_TEST_CASE("SAD", "", test_scf::float_types) {
    using float_type     = TestType;
    using allocator_type = tensorwrapper::allocator::Eigen<float_type>;
    using rank2_il_type  = typename allocator_type::rank2_il;

    auto mm  = test_scf::load_modules<float_type>();
    auto aos = test_scf::h2_aos();
    auto H   = test_scf::h2_hamiltonian();
    auto rt  = mm.get_runtime();

    auto mod = mm.at("SAD guess");
    auto psi          = mod.template run_as<pt>(H, aos);
    const auto& evals = psi.orbitals().diagonalized_matrix();

    occ_index occs{0};
    allocator_type alloc(rt);
    shape_type shape_corr{2};
    auto pbuffer = alloc.construct({-0.498376, 0.594858});
    tensor corr(shape_corr, std::move(pbuffer));

    REQUIRE(psi.orbital_indices() == occs);
    REQUIRE(psi.orbitals().from_space() == aos);
    REQUIRE(approximately_equal(corr, evals, 1E-6));
}