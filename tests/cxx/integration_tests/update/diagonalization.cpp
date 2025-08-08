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

TEMPLATE_LIST_TEST_CASE("Diagaonalization", "", test_scf::float_types) {
    using float_type = TestType;
    auto mm          = test_scf::load_modules<float_type>();
    auto& mod        = mm.at("Diagonalization Fock update");

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

        const auto& psi = mod.template run_as<pt>(f_e, core_guess);
        REQUIRE(psi.orbital_indices() == occs);
        REQUIRE(psi.orbitals().from_space() == aos);

        // Check orbital energies
        const auto& evals = psi.orbitals().diagonalized_matrix();
        tensorwrapper::allocator::Eigen<float_type> alloc(mm.get_runtime());
        auto corr_buffer = alloc.construct({-1.25330893, -0.47506974});
        tensorwrapper::shape::Smooth corr_shape{2};
        tensorwrapper::Tensor corr(corr_shape, std::move(corr_buffer));

        using tensorwrapper::operations::approximately_equal;
        REQUIRE(approximately_equal(evals, corr, 1E-6));
    }
}
