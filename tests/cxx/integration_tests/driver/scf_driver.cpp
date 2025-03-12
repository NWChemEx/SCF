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

using pt = simde::AOEnergy;

TEMPLATE_LIST_TEST_CASE("SCFDriver", "", test_scf::float_types) {
    using float_type = TestType;
    auto mm          = test_scf::load_modules<float_type>();
    auto h2          = test_scf::make_h2<simde::type::chemical_system>();
    auto aos         = test_scf::h2_aos().ao_basis_set();

    tensorwrapper::allocator::Eigen<float_type> alloc(mm.get_runtime());
    tensorwrapper::shape::Smooth shape_corr{};
    auto pcorr = alloc.allocate(tensorwrapper::layout::Physical(shape_corr));
    using tensorwrapper::operations::approximately_equal;
    pcorr->at() = -1.1167592336;
    tensorwrapper::Tensor corr(shape_corr, std::move(pcorr));

    const auto e = mm.template run_as<pt>("SCF Driver", aos, h2);
    REQUIRE(approximately_equal(corr, e, 1E-6));

    SECTION("H2 Dimer") {
        simde::type::nucleus h0("H", 1ul, 1836.15, 0.0, 0.0, 0.0);
        simde::type::nucleus h1("H", 1ul, 1836.15, 0.0, 0.0, 1.39839);
        simde::type::nucleus h2("H", 1ul, 1836.15, 0.0, 0.0, 4.39839);
        simde::type::nucleus h3("H", 1ul, 1836.15, 0.0, 0.0, 5.79678);
        simde::type::nuclei h2_dimer_nuclei{h0, h1, h2, h3};
        auto ao_bs = test_scf::h_basis(h2_dimer_nuclei);
        simde::type::molecule h2_dimer_mol(0, 1, h2_dimer_nuclei);
        simde::type::chemical_system h2_dimer_sys(h2_dimer_mol);
        const auto e =
          mm.template run_as<pt>("SCF Driver", ao_bs, h2_dimer_sys);
        alloc.rebind(corr.buffer()).at() = -2.2260535919670001;
        REQUIRE(approximately_equal(corr, e, 1E-6));
    }
}